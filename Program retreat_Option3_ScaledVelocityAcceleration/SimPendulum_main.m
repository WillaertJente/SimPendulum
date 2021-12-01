%% SimPendulum_main

%% Input
% Change here subject and trial
info.subj   = 'TD5';           % Subject name
info.trial  = 1;               % Trial number
info.option = 'Opt3_Scaled';              % Name to save results
info.wq     = 1;               % weight on q error
info.wqd    = 0.5;             % weight on qd error

%% Import subject parameters and experimental data

% Path info -  Path to model and experimental data
pathmain = pwd;
[pathTemp,~,~] = fileparts(pathmain);
[pathRepo,~,~] = fileparts(pathTemp);
info.path      = [pathTemp '\Implicit\Muscle\Experimental data\' info.subj '\'];

% Import subject parameters
params_subject = ImportSubjectParameters(info);    % Input parameters (mtot,lc, l, age, m, RG, SE, Nmr, z)

% Import experimental data
bool_plot = 0;
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
x          = opti.variable(1,N);
xd_scaled  = opti.variable(1,N);  % x * 5
xdd_scaled = opti.variable(1,N);  % x * 25
lMtilda    = opti.variable(1,N);

% (Slack) controls
vMtilda      = opti.variable(1,N);
lM_projected = opti.variable(1,N);

% Parameters that will be optimized
a            = opti.variable(1);
kFpe         = opti.variable(1);
B            = opti.variable(1); 

% Bounds
[Ub Lb] = SelectBounds();
opti.subject_to(-pi     <  x    < Ub.x);
opti.subject_to(Lb.a    <  a    < Ub.a);
opti.subject_to(Lb.kFpe <  kFpe < Ub.kFpe);    % nominal value = 0.1
opti.subject_to(Lb.B    < B     < Ub.B); 
opti.subject_to(Lb.lM_projected  <  lM_projected);
opti.subject_to(Lb.vMtilda  <  vMtilda < Ub.vMtilda);
opti.subject_to(Lb.lMtilda  <  lMtilda < Ub.lMtilda);
% opti.subject_to(min(data_exp.qdspline)-2  < xd  < max(data_exp.qdspline)+2);

% Constraints on initial states
opti.subject_to(x(1)            == data_exp.x0(1));
opti.subject_to(xd_scaled(1)    == 0);

% Initial guess
opti.set_initial(x, data_exp.qspline);
opti.set_initial(xd_scaled*5,data_exp.qdspline);
opti.set_initial(a, 0.01);
opti.set_initial(kFpe,0.1);
opti.set_initial(B,0.1); 
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
[error] = CalculateMusculoSkeletalDynamics(x,xd_scaled*5,xdd_scaled*25, lMtilda, lM_projected, kFpe,vMtilda,a, data_exp, coeff_LMT_ma, params_OS, shift, B); 

% Constraints - trapezoidal integration 
dt = dt_spline; 
opti.subject_to((xd_scaled(1:N-1)*5 + xd_scaled(2:N)*5)*dt/2 + x(1:N-1) == x(2:N));
opti.subject_to((xdd_scaled(2:N)*25+xdd_scaled(1:N-1)*25)*dt/2 +xd_scaled(1:N-1)*5 == xd_scaled(2:N)*5);
opti.subject_to((dlMdt(2:N)+dlMdt(1:N-1))*dt/2 + lMtilda(1:N-1) == lMtilda(2:N)); % VmTilde met factor 10 (MRS)
opti.subject_to(error == 0);

% Objective function
J = DefineObjectiveFunction(x,xd_scaled*5,data_exp, info); 
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
R.info   = info; 
R.subject= params_subject;
R.OS     = params_OS;

% States
R.x       = sol.value(x); R.xd = sol.value(xd_scaled); R.xdd = sol.value(xdd_scaled); 
R.lMtilda = sol.value(lMtilda); 

% Controls
R.vMtilda     = sol.value(vMtilda);
R.lMprojected = sol.value(lM_projected); 
R.dlMdt       = sol.value(dlMdt); 

% Parameters
R.a       = sol.value(a); 
R.kFpe    = sol.value(kFpe);
R.B       = sol.value(B); 

% Objective function
R.wq      = info.wq;
R.wqd     = info.wqd;
R.J       = sol.value(J);

% Additional information
R.error   = sol.value(error); 
R.coeff   = coeff_LMT_ma; 
R.initGuess= InitGuess; 
R.shift   = shift; 

% Bounds
R.bounds.Ub = Ub; 
R.bounds.Lb = Lb; 

% Calculated Parameters
[R.C] = RecalculateOutcomes(R); 

% Write results 
save([pathTemp,'/Results/',info.subj,'_T',num2str(info.trial),'_',info.option,'.mat'],'R')

%% Forward version
[q_forward,qd_forward,lMtilda_forward] = forwardSim(R.x(1),R.xd(1),R.lMtilda(:,1),R.kFpe,R.a ,R.exp, coeff_LMT_ma, params_OS, shift, dt, N);

%% Plot
% Results
q_forward = []; 
h = PlotResults(R, q_forward, info);
saveas(h,[pathTemp,'/Results/Figures/', info.subj, '_T',num2str(info.trial),info.option,'_1.Results.fig']);

% Params
g = PlotParams(R,info);
saveas(g,[pathTemp,'/Results/Figures/', info.subj, '_T',num2str(info.trial),info.option,'_2.Params.fig']);

% Muscle geometry
f = PlotMuscleGeometry(R, info);
saveas(f,[pathTemp,'/Results/Figures/', info.subj, '_T',num2str(info.trial),info.option,'_3.MuscleGeometry.fig']);

% Torques
p = PlotTorques(R); 
saveas(f,[pathTemp,'/Results/Figures/', info.subj, '_T',num2str(info.trial),info.option,'_4.Torques.fig']);
