%% SimPendulum_main
%% Input
% Change here subject and trial
info.subj   = 'TD5';           % Subject name
info.trial  = 2;               % Trial number
info.option = 'Opt5_SRS_FixedTime';              % Name to save results
info.wq     = 1;               % weight on q error
info.wqd    = 0.5;             % weight on qd error
info.kSRS   = 0; 

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
N_1     = data_exp.N_1; 
% States
x            = opti.variable(1,N);
xd           = opti.variable(1,N); 
xdd          = opti.variable(1,N); 
lMtilda      = opti.variable(2,N);  
Fsrs         = opti.variable(1,N-N_1); 

% (Slack) controls
vMtilda      = opti.variable(2,N);
dFsrsdt      = opti.variable(1,N-N_1); 
lM_projected = opti.variable(2,N);

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
opti.subject_to(Lb.vMtilda  <  vMtilda  < Ub.vMtilda);
opti.subject_to(Lb.lMtilda  <  lMtilda  < Ub.lMtilda);

% % Calculate initial value of Fsrs 
[Fsrs_f1_cal] = CalculateInitialValueFsrs(lMtilda, a_ext, info, params_OS, data_exp); 

% Constraints on initial states
opti.subject_to(x(1)     == data_exp.x0(1));
opti.subject_to(xd(1)    == 0);
opti.subject_to(Fsrs(1)  == Fsrs_f1_cal);  

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

% Skeletal dynamics fase 1 
[error_f1] = CalculateMusculoSkeletalDynamics_F1(x(1:N_1),xd(1:N_1),xdd(1:N_1), lMtilda(:,1:N_1), lM_projected(:,1:N_1),kFpe,vMtilda(:,1:N_1), a_ext, a_flex, data_exp, coeff_LMT_ma, params_OS, shift, B, info); 

% Skeletal dynamics fase 2
[error_f2] = CalculateMusculoSkeletalDynamics_F2(x(N_1+1:end),xd(N_1+1:end),xdd(N_1+1:end),lMtilda(:,N_1+1:end), lM_projected(:,N_1+1:end), kFpe, vMtilda(:,N_1+1:end), ...
             a_ext, a_flex, data_exp, coeff_LMT_ma, params_OS, shift, B, info, Fsrs, dFsrsdt); 

% Constraints - trapezoidal integration 
dt = dt_spline; 
lMtilda_ext = lMtilda(1,:); 
lMtilda_flex= lMtilda(2,:); 
opti.subject_to((xd(1:N-1) + xd(2:N))*dt/2 + x(1:N-1) == x(2:N));
opti.subject_to((xdd(2:N)+xdd(1:N-1))*dt/2 +xd(1:N-1) == xd(2:N));
opti.subject_to((dlMdt_ext(2:N)+dlMdt_ext(1:N-1))*dt/2 + lMtilda_ext(1:N-1) == lMtilda_ext(2:N));
opti.subject_to((dlMdt_flex(2:N)+dlMdt_flex(1:N-1))*dt/2 + lMtilda_flex(1:N-1) == lMtilda_flex(2:N));
opti.subject_to((dFsrsdt(2:N-N_1)+dFsrsdt(1:N-N_1-1))*dt/2 + Fsrs(1:N-N_1-1) == Fsrs(2:N-N_1));     % Constraint on variable, in functie error dat constraint = variable
opti.subject_to(error_f1 == 0);
opti.subject_to(error_f2 == 0); 
% opti.subject_to(Fsrs(1) == Fsrs_f1_cal) % Fsrs start fase 2 moet gelijk zijn aan Fsrs einde fase 1 

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
R.error_f1  = sol.value(error_f1);
R.error_f2  = sol.value(error_f2);
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
%[q_forward,qd_forward,lMtilda_forward] = forwardSim(R.x(1),R.xd(1),R.lMtilda(:,1),R.kFpe,R.a ,R.exp, coeff_LMT_ma, params_OS, shift, dt, N, R.B);
[q_forward_F1,qd_forward_F1,lMtilda_forward_F1] = forwardSim_F1(R.x(1),R.xd(1),R.lMtilda(:,1),R.kFpe,R.a ,R.exp, coeff_LMT_ma, params_OS, shift, dt, data_exp.N_1-1, R.B, info);
[q_forward_F2,qd_forward_F2,lMtilda_forward_F2] = forwardSim_F2(R.x(N_1+1),R.xd(N_1+1),R.lMtilda(:,1),R.kFpe,R.a ,R.exp, coeff_LMT_ma, params_OS, shift, dt, data_exp.N_1+1:end, R.B, info);


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
