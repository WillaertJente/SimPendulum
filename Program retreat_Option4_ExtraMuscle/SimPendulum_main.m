%% SimPendulum_main
%% Input
% Change here subject and trial
info.subj   = 'TD5';           % Subject name
info.trial  = 1;               % Trial number
info.option = 'Opt4_AddMuscle_AddVmtoJ';              % Name to save results
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
x            = opti.variable(1,N);
xd           = opti.variable(1,N); 
xdd          = opti.variable(1,N); 
lMtilda      = opti.variable(2,N);  % RF
%lMtilda_flex = opti.variable(1,N);  % BFlh

% (Slack) controls
vMtilda      = opti.variable(2,N);
%vMtilda_flex     = opti.variable(1,N);
lM_projected = opti.variable(2,N);
%lM_projected_flex= opti.variable(1,N); 

% Parameters that will be optimized
a_ext         = opti.variable(1);
a_flex        = opti.variable(1); 
kFpe          = opti.variable(1);
B             = opti.variable(1); 

% Bounds
[Ub, Lb] = SelectBounds();
opti.subject_to(Lb.x   <  x     < Ub.x);
opti.subject_to(Lb.a   < a_ext  < Ub.a);
opti.subject_to(Lb.a   < a_flex < Ub.a);
opti.subject_to(Lb.kFpe <  kFpe < Ub.kFpe);    % nominal value = 0.1
opti.subject_to(Lb.B    < B     < Ub.B); 
opti.subject_to(Lb.lM_projected  <  lM_projected(1,:));
opti.subject_to(Lb.lM_projected  <  lM_projected(2,:));
%opti.subject_to(Lb.lM_projected  <  lM_projected_flex);
opti.subject_to(Lb.vMtilda  <  vMtilda  < Ub.vMtilda);
%opti.subject_to(Lb.vMtilda  <  vMtilda_flex < Ub.vMtilda);
opti.subject_to(Lb.lMtilda  <  lMtilda  < Ub.lMtilda);
%opti.subject_to(Lb.lMtilda  <  lMtilda_flex < Ub.lMtilda);
% opti.subject_to(min(data_exp.qdspline)-2  < xd  < max(data_exp.qdspline)+2);

% Constraints on initial states
opti.subject_to(x(1)     == data_exp.x0(1));
opti.subject_to(xd(1)    == 0);

% Initial guess
opti.set_initial(x, data_exp.qspline);
opti.set_initial(xd,data_exp.qdspline);
opti.set_initial(a_ext, 0.01);
opti.set_initial(a_flex, 0.01);
opti.set_initial(kFpe,0.1);
opti.set_initial(B,0.1); 
opti.set_initial(lM_projected, InitGuess.lM_projected');
opti.set_initial(lMtilda, InitGuess.lMtilda');     
opti.set_initial(vMtilda, InitGuess.vM);



%% Define problem (muscle model)
% Calculate shift
kT     = 35;
shift  = getshift(kT);

% Calculate dlMdt
[dlMdt]   = CalculateDLMDT(vMtilda, params_OS); 
dlMdt_ext = dlMdt(1,:); dlMdt_flex = dlMdt(2,:); 

% Skeletal dynamics 
[error] = CalculateMusculoSkeletalDynamics(x,xd,xdd, lMtilda, lM_projected,kFpe,vMtilda, a_ext, a_flex, data_exp, coeff_LMT_ma, params_OS, shift, B); 

% Constraints - trapezoidal integration 
dt = dt_spline; 
lMtilda_ext = lMtilda(1,:); 
lMtilda_flex= lMtilda(2,:); 
opti.subject_to((xd(1:N-1) + xd(2:N))*dt/2 + x(1:N-1) == x(2:N));
opti.subject_to((xdd(2:N)+xdd(1:N-1))*dt/2 +xd(1:N-1) == xd(2:N));
opti.subject_to((dlMdt_ext(2:N)+dlMdt_ext(1:N-1))*dt/2 + lMtilda_ext(1:N-1) == lMtilda_ext(2:N));
opti.subject_to((dlMdt_flex(2:N)+dlMdt_flex(1:N-1))*dt/2 + lMtilda_flex(1:N-1) == lMtilda_flex(2:N));
opti.subject_to(error == 0);

% Objective function
J = DefineObjectiveFunction(x,xd,data_exp, info, vMtilda); 
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
R.x           = sol.value(x); R.xd = sol.value(xd); R.xdd = sol.value(xdd); 
R.lMtilda     = sol.value(lMtilda); 

% Controls
R.vMtilda = sol.value(vMtilda);
R.lMprojected = sol.value(lM_projected);
R.dlMdt            = sol.value(dlMdt); 

% Parameters
R.a_ext   = sol.value(a_ext);
R.a_flex  = sol.value(a_flex);
R.a       = [R.a_ext R.a_flex]; 
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
R.a(1) = 0.005; 
[q_forward,qd_forward,lMtilda_forward] = forwardSim(R.x(1),R.xd(1),R.lMtilda(:,1),R.kFpe,R.a ,R.exp, coeff_LMT_ma, params_OS, shift, dt, N, R.B);

%% Plot
% Results
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
