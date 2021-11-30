%% SimPendulum_main

%% Input
% Change here subject and trial
info.subj   = 'TD5';           % Subject name
info.trial  = 1;               % Trial number
info.option = '';              % Name to save results
info.wq     = 1;               % weight on q error
info.wqd    = 0;             % weight on qd error

%% Import subject parameters and experimental data

% Path info -  Path to model and experimental data
pathmain = pwd;
[pathTemp,~,~] = fileparts(pathmain);
[pathRepo,~,~] = fileparts(pathTemp);
info.path      = [pathTemp '\Implicit\Muscle\Experimental data\' info.subj '\'];

% Import subject parameters
params_subject = ImportSubjectParameters(info);    % Input parameters (mtot,lc, l, age, m, RG, SE, Nmr, z)

% Import experimental data
bool_plot = 1;
dt_spline = 0.005;
data_exp  = ImportExperimentalData(info, bool_plot, params_subject, dt_spline);

% Define phases of pendulum (initial state, end of first swing)
bool_plot = 1;
[data_exp.x0, data_exp.N_1] = PendulumPhases(data_exp, bool_plot);

%% Import OpenSim model parameters
% OpenSim
addpath('MuscleModel');
muscles = {'rect_fem_','bifemlh_'};
[params_OS] = ReadOpenSimParams(info, params_subject, muscles);

%% Calculate LMT en Ma coefficients
bool_plot = 1;
[coeff_LMT_ma] = DefineLMTCoefficients(params_subject, info, muscles, bool_plot);

%% Create initial guess
bool_guess  = 1; % to create trial specific initial guess of lMtilda
[InitGuess] = CreateInitialGuess(bool_guess,params_OS, data_exp, muscles, coeff_LMT_ma);


%% Define states, controls, bounds and initial guess
import casadi.*;        % Import casadi libraries
opti = casadi.Opti();   % Initialise opti structure

N       = data_exp.Nspline;
% States
x       = opti.variable(1,N);
xd      = opti.variable(1,N);
xdd     = opti.variable(1,N); 
lMtilda = opti.variable(1,N);

% (Slack) controls
vMtilda      = opti.variable(1,N);
lM_projected = opti.variable(1,N);

% Parameters that will be optimized
a            = opti.variable(1);
kFpe         = opti.variable(1);

% Bounds
opti.subject_to(-pi   <  x  < 10*pi/180);
% opti.subject_to(min(data_exp.qdspline)-2  < xd  < max(data_exp.qdspline)+2);
opti.subject_to(0.001 <  a   < 1);
opti.subject_to(1e-4  <  lM_projected);
opti.subject_to(-10   <  vMtilda < 10);
opti.subject_to(0.2   <  lMtilda < 1.5);
opti.subject_to(0     <  kFpe    < 0.2);    % nominal value = 0.1
% opti.subject_to(0.001==  a   );

% Constraints on initial states
opti.subject_to(x(1)     == data_exp.x0(1));
opti.subject_to(xd(1)    == 0);

% Initial guess
opti.set_initial(x, data_exp.qspline);
opti.set_initial(xd,data_exp.qdspline);
opti.set_initial(a, 0.01);
opti.set_initial(kFpe,0.1);
opti.set_initial(lM_projected, InitGuess.lM_projected(:,1));
opti.set_initial(lMtilda, InitGuess.lMtilda(:,1));         % LmTildeGuess
opti.set_initial(vMtilda, InitGuess.vM(1,:));

%% Define problem (muscle model)
% Calculate shift
kT     = 35;
shift  = getshift(kT);

% Calculate dlMdt
[dlMdt] = CalculateDLMDT(vMtilda, params_OS); 

% Skeletal dynamics 
[error] = CalculateMusculoSkeletalDynamics(x,xd,xdd, lMtilda, lM_projected, kFpe,vMtilda,a, data_exp, coeff_LMT_ma, params_OS, shift); 

% Constraints - forward euler 
dt = dt_spline; 
opti.subject_to((xd(1:N-1) + xd(2:N))*dt/2 + x(1:N-1) == x(2:N));
opti.subject_to((xdd(2:N)+xdd(1:N-1))*dt/2 +xd(1:N-1) == xd(2:N));
opti.subject_to((dlMdt(2:N)+dlMdt(1:N-1))*dt/2 + lMtilda(1:N-1) == lMtilda(2:N)); % VmTilde met factor 10 (MRS)
opti.subject_to(error == 0);

% Objective function
J = DefineObjectiveFunction(x,xd,data_exp, info); 
opti.minimize(J); 
    
%% Solve problem
% options for IPOPT
options.ipopt.tol = 1*10^(-6);          
options.ipopt.linear_solver = 'mumps';
          
% Solve the OCP
opti.solver('ipopt',options);
sol = opti.solve();  
    
%% Results 
% Experimental data 
R.exp    = data_exp;

% States
R.x       = sol.value(x); R.xd = sol.value(xd); R.xdd = sol.value(xdd); 
R.lMtilda = sol.value(lMtilda); 

% Controls
R.vMtilda = sol.value(vMtilda);
R.lMprojected = sol.value(lM_projected); 

% Parameters
R.a       = sol.value(a); 
R.kFpe    = sol.value(kFpe);

% Objective function
R.wq      = info.wq;
R.wqd     = info.wqd;
R.J       = sol.value(J);

% Skelet dynamica
% R.skelet.error_force = sol.value(error_force);
% R.skelet.dlMdt       = sol.value(dlMdt);
% R.skelet.error_skelet= sol.value(error_skelet);
% R.skelet.lM          = sol.value(lM); 
% R.skelet.w           = sol.value(w);

% Forces
% R.FT   = sol.value(FT);
% R.Fpe  = sol.value(Fpe);
% R.Fce  = sol.value(Fce); 

% Muscle tendon length and moment arms
% R.lMT  = sol.value(lMT); 
% R.MA   = sol.value(MA); 
% 
% R.T = R.FT.*R.MA;
%% Plot
figure()
plot(data_exp.qspline*180/pi,'k','LineWidth',1.5); hold on
plot(R.x*180/pi,'r','LineWidth',1.5); hold on; 
box off; legend({'Exp','Sim'}); ylabel('Knee Angle ({\circ})'); xlabel('Frames')

figure()
subplot(421)
plot(R.exp.qspline*180/pi,'k','LineWidth',1.5); 
subplot(422)
plot(R.x*180/pi,'k','LineWidth',1.5); hold on;
title('Solution')
subplot(423)
plot(R.lMtilda,'k','LineWidth',1.5); hold on
title('lMtilda');
subplot(426)
plot(R.vMtilda,'k','LineWidth',1.5); hold on
title('vMtilda')
subplot(424)
plot(R.lMprojected,'k','LineWidth',1.5); hold on
title('lMprojected')
subplot(425)
plot(R.skelet.lM,'k','LineWidth',1.5); hold on
title('lM')
subplot(427)
plot(R.skelet.dlMdt,'k','LineWidth',1.5); hold on
title('dLmdt')
subplot(428)
plot(R.skelet.w,'k','LineWidth',1.5); hold on

figure()
subplot(311)
plot(R.FT,'k','LineWidth',1.5); hold on
title('FT')
subplot(312)
plot(R.Fce,'k','LineWidth',1.5); hold on
title('Fce')
subplot(313)
plot(R.Fpe,'k','LineWidth',1.5); hold on
title('FPe')

figure()
subplot(121)
bar(1,R.kFpe); hold on
title('Passive stiffness')
subplot(122)
bar(1,R.a); hold on
title('activation')

figure()
subplot(211)
plot(R.lMT(1,:),'k','LineWidth',1.5); 
title('lMT')
subplot(212)
plot(R.MA(1,:),'k','LineWidth',1.5); 
title('MA')