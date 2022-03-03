%% SimPendulum_main
%% Input
% Change here subject and trial
info.subj   = 'CP10';           % Subject name
info.option = 'Opt8_MT';              % Name to save results
info.wq     = 1;               % weight on q error
info.wqd    = 0.5;             % weight on qd error
info.kSRS   = 280;
info.tau    = 0.080;
n_tr        = 2;
% B = 0.0467;

%% Import subject parameters
% Path info -  Path to model and experimental data
pathmain = pwd;
[pathTemp,~,~] = fileparts(pathmain);
[pathRepo,~,~] = fileparts(pathTemp);
info.path      = [pathTemp '\Implicit\Muscle\Experimental data\New/' info.subj '\'];

% Import subject parameters
params_subject = ImportSubjectParameters(info);    % Input parameters (mtot,lc, l, age, m, RG, SE, Nmr, z)

%% Import experimental data
% Trial 1
info.trial  = [12];
% Import experimental data
bool_plot = 0;
dt_spline = 0.005;
data_exp  = ImportExperimentalData(info, bool_plot, params_subject, dt_spline);

% Define phases of pendulum (initial state, end of first swing)
bool_plot = 1;
[data_exp.x0, data_exp.N_1] = PendulumPhases(data_exp, bool_plot);

% Part 1 data
data_p1 = data_exp;

% Trial 2
info.trial  = [13];
% Import experimental data
bool_plot = 0;
dt_spline = 0.005;
data_exp  = ImportExperimentalData(info, bool_plot, params_subject, dt_spline);

% Define phases of pendulum (initial state, end of first swing)
bool_plot = 1;
[data_exp.x0, data_exp.N_1] = PendulumPhases(data_exp, bool_plot);

% Part 2 data
data_p2 = data_exp;

% Combine data
clear data_exp;
[data_exp] = CreateTotalData(data_p1, data_p2);

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

%% Define start and stop for total vector
start    = [1 data_exp.Nspline(1)+1];
stop     = [data_exp.Nspline(1) data_exp.Nspline(1)+data_exp.Nspline(2)];
start_N1 = [1 data_exp.Nspline(1)-data_exp.N_1(1)+1];
stop_N1  = [data_exp.Nspline(1)-data_exp.N_1(1) data_exp.Nspline(1)-data_exp.N_1(1) + data_exp.Nspline(2)-data_exp.N_1(2)];

%%
import casadi.*;        % Import casadi libraries
opti = casadi.Opti();   % Initialise opti structure

%% Define vector for casadi variables with length of all trials together
N_tot            = data_exp.Nspline(1) + data_exp.Nspline(2);
N_phase          = (data_exp.Nspline(1)-data_exp.N_1(1)) + (data_exp.Nspline(2)-data_exp.N_1(2));
% States
x_tot            = opti.variable(1,N_tot);
xd_tot           = opti.variable(1,N_tot);
xdd_tot          = opti.variable(1,N_tot);
lMtilda_tot      = opti.variable(2,N_tot);
Fsrs_tot         = opti.variable(1,N_phase);
Fsrs_del_tot     = opti.variable(1,N_tot);

% (Slack) controls
vMtilda_tot      = opti.variable(2,N_tot);
dFsrsdt_tot      = opti.variable(1,N_phase);
dFsrs_deldt_tot  = opti.variable(1,N_tot);
lM_projected_tot = opti.variable(2,N_tot);
dt1_tot          = opti.variable(2);

% Parameters that will be optimized
a_ext_tot        = opti.variable(2);
a_flex_tot       = opti.variable(2);
kFpe_tot         = opti.variable(1);
B_tot            = opti.variable(1);
kR_tot           = opti.variable(2);

%% Loop over all trials
J = 0;
for i = 1:2
    N_1 = data_exp.N_1(i);
    N   = data_exp.Nspline(i); 
    
    % States
    x            = x_tot(1,start(i):stop(i));
    xd           = xd_tot(1,start(i):stop(i));
    xdd          = xdd_tot(1,start(i):stop(i));
    lMtilda      = lMtilda_tot(:,start(i):stop(i));
    Fsrs         = Fsrs_tot(1,start_N1(i):stop_N1(i));
    Fsrs_del     = Fsrs_del_tot(1,start(i):stop(i));
    
    % (Slack) controls
    vMtilda      = vMtilda_tot(:,start(i):stop(i));
    dFsrsdt      = dFsrsdt_tot(1,start_N1(i):stop_N1(i));
    dFsrs_deldt  = dFsrs_deldt_tot(1,start(i):stop(i));
    lM_projected = lM_projected_tot(:,start(i):stop(i));
    dt1          = dt1_tot(i);
    
    % Parameters that will be optimized
    a_ext         = a_ext_tot(i);
    a_flex        = a_flex_tot(i);
    kFpe          = kFpe_tot(1);
    B             = B_tot(1);
    kR            = kR_tot(i);
    
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
    opti.subject_to(Lb.dt1 < dt1 < Ub.dt1);
    opti.subject_to(Lb.kR  < kR  < Ub.kR);
    
    % Calculate initial value of Fsrs
    [Fsrs_f1_cal] = CalculateInitialValueFsrs(lMtilda, a_ext, info, params_OS, N_1);
    
    % Constraints on initial states
    opti.subject_to(x(1)     == data_exp.x0(i,1));
    opti.subject_to(xd(1)    == 0);
    opti.subject_to(Fsrs(1)  == Fsrs_f1_cal);
    
    % Initial guess
    opti.set_initial(x, data_exp.qspline(start(i):stop(i)));
    opti.set_initial(xd,data_exp.qdspline(start(i):stop(i)));
    opti.set_initial(a_ext, 0.01);
    opti.set_initial(a_flex, 0.01);
    opti.set_initial(kFpe,0.1);
    opti.set_initial(B,0.1);
    opti.set_initial(lM_projected, InitGuess.lM_projected(start(i):stop(i),:)');
    opti.set_initial(lMtilda, InitGuess.lMtilda(start(i):stop(i),:)');
    opti.set_initial(vMtilda, InitGuess.vM(:,start(i):stop(i)));
    opti.set_initial(dt1,0.005);
    opti.set_initial(kR, 0.01);
    
    %% Define problem (muscle model)
    % Calculate shift
    kT     = 35;
    shift  = getshift(kT);
    
    % Calculate dlMdt
    [dlMdt]   = CalculateDLMDT(vMtilda, params_OS);
    dlMdt_ext = dlMdt(1,:); dlMdt_flex = dlMdt(2,:);
    
    % Time phase 1
    tF1 = N_1*dt1;
    dt2 = 0.005;
    
    opti.subject_to(tF1 > 0.2);
    opti.subject_to(xd(1:N_1) < 1e-4);
    opti.subject_to(1e-4 < xd(N_1+1)); % of kleiner als 1 er nog bij?
    
    % Skeletal dynamics fase 1
    lMtilda_init = lMtilda(1,1);
    [error_f1] = CalculateMusculoSkeletalDynamics_F1_Refl(x(1:N_1),xd(1:N_1),xdd(1:N_1), lMtilda(:,1:N_1), lMtilda_init, lM_projected(:,1:N_1),kFpe,vMtilda(:,1:N_1), a_ext, a_flex, data_exp, coeff_LMT_ma, params_OS, shift, B, info, Fsrs_del(1:N_1), dFsrs_deldt(1:N_1), kR);
    
    % Skeletal dynamics fase 2
    [error_f2] = CalculateMusculoSkeletalDynamics_F2_Refl(x(N_1+1:end),xd(N_1+1:end),xdd(N_1+1:end),lMtilda(:,N_1+1:end), lM_projected(:,N_1+1:end), kFpe, vMtilda(:,N_1+1:end), ...
        a_ext, a_flex, data_exp, coeff_LMT_ma, params_OS, shift, B, info, Fsrs, dFsrsdt,Fsrs_del(N_1+1:end), dFsrs_deldt(N_1+1:end), kR);
    
    % Constraints - trapezoidal integration
    % dt = dt_spline;
    lMtilda_ext = lMtilda(1,:);
    lMtilda_flex= lMtilda(2,:);
    opti.subject_to((xd(1:N_1-1) + xd(2:N_1))*dt1/2 + x(1:N_1-1) == x(2:N_1));
    opti.subject_to((xd(N_1:N-1) + xd(N_1+1:N))*dt2/2 + x(N_1:N-1) == x(N_1+1:N));
    opti.subject_to((xdd(2:N_1)+xdd(1:N_1-1))*dt1/2 +xd(1:N_1-1) == xd(2:N_1));
    opti.subject_to((xdd(N_1+1:N)+xdd(N_1:N-1))*dt2/2 +xd(N_1:N-1) == xd(1+N_1:N));
    opti.subject_to((dlMdt_ext(2:N_1)+dlMdt_ext(1:N_1-1))*dt1/2 + lMtilda_ext(1:N_1-1) == lMtilda_ext(2:N_1));
    opti.subject_to((dlMdt_ext(1+N_1:N)+dlMdt_ext(N_1:N-1))*dt2/2 + lMtilda_ext(N_1:N-1) == lMtilda_ext(1+N_1:N));
    opti.subject_to((dlMdt_flex(2:N_1)+dlMdt_flex(1:N_1-1))*dt1/2 + lMtilda_flex(1:N_1-1) == lMtilda_flex(2:N_1));
    opti.subject_to((dlMdt_flex(1+N_1:N)+dlMdt_flex(N_1:N-1))*dt2/2 + lMtilda_flex(N_1:N-1) == lMtilda_flex(1+N_1:N));
    opti.subject_to((dFsrsdt(2:N-N_1)+dFsrsdt(1:N-N_1-1))*dt2/2 + Fsrs(1:N-N_1-1) == Fsrs(2:N-N_1));     % Constraint on variable, in functie error dat constraint = variable
    opti.subject_to((dFsrs_deldt(1:N_1-1) + dFsrs_deldt(2:N_1))*dt1/2 + Fsrs_del(1:N_1-1) == Fsrs_del(2:N_1));
    opti.subject_to((dFsrs_deldt(N_1:N-1) + dFsrs_deldt(N_1+1:N))*dt2/2 + Fsrs_del(N_1:N-1) == Fsrs_del(N_1+1:N));
    opti.subject_to(error_f1 == 0);
    opti.subject_to(error_f2 == 0);
    
    % Objective function
    qspline  = data_exp.qspline(start(i):stop(i));
    qdspline = data_exp.qdspline(start(i):stop(i));
    
    cost = DefineObjectiveFunction(x, xd, qspline, qdspline, info, vMtilda);
    J    = J + cost;
end

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
R.x           = sol.value(x_tot); R.xd = sol.value(xd_tot); R.xdd = sol.value(xdd_tot);
R.lMtilda     = sol.value(lMtilda_tot);

% Controls
R.vMtilda     = sol.value(vMtilda_tot);
R.lMprojected = sol.value(lM_projected_tot);
R.Fsrs        = sol.value(Fsrs_tot);
R.Fsrs_del    = sol.value(Fsrs_del_tot);
R.dFsrs_deldt = sol.value(dFsrs_deldt_tot);

% Parameters
R.a_ext   = sol.value(a_ext_tot);
R.a_flex  = sol.value(a_flex_tot);
R.kFpe    = sol.value(kFpe_tot);
R.B       = sol.value(B_tot);
R.dt1     = sol.value(dt1_tot);
R.dt2     = dt2;
R.kR      = sol.value(kR_tot);

% Objective function
R.wq      = info.wq;
R.wqd     = info.wqd;
R.J       = sol.value(J);

% Additional information
R.error_f1  = sol.value(error_f1);
R.error_f2  = sol.value(error_f2);
R.coeff     = coeff_LMT_ma;
R.initGuess = InitGuess;
R.shift     = shift;

% Bounds
R.bounds.Ub = Ub;
R.bounds.Lb = Lb;

% Calculated Parameters
%[R.C] = RecalculateOutcomes2(R, info, N_1);

% Write results
save([pathTemp,'/Results/',info.subj,'_T',num2str(info.trial),'_',info.option,'.mat'],'R')

% %% Forward version
% %[q_forward,qd_forward,lMtilda_forward] = forwardSim(R.x(1),R.xd(1),R.lMtilda(:,1),R.kFpe,R.a ,R.exp, coeff_LMT_ma, params_OS, shift, dt, N, R.B);
% % F1
% [q_forward_F1,qd_forward_F1,lMtilda_forward_F1, Fsrs_del_forward_F1] = forwardSim_F1_Refl(R.x(1),R.xd(1),R.lMtilda(:,1),R.kFpe,R.a ,R.exp, coeff_LMT_ma, params_OS, shift, R.dt1, data_exp.N_1-1, R.B, info, R.Fsrs_del(1), R.kR);
% [Fsrs_F1] =  CalculateFsrsF1Forward(lMtilda_forward_F1(1,:),R.lMtilda(1,1), R.a(1), info.kSRS, params_OS);
% 
% % F2
% %[Fsrs] = CalculateFsrsForward(lMtilda_forward_F1(1,:),R.lMtilda(1,1), R.a(1), info.kSRS, params_OS);
% [q_forward_F2,qd_forward_F2,lMtilda_forward_F2, Fsrs_forward_F2, Fsrs_del_forward_F2] = forwardSim_F2_Refl(q_forward_F1(end),qd_forward_F1(end),lMtilda_forward_F1(:,end),R.kFpe,R.a ,R.exp, coeff_LMT_ma, params_OS, shift, dt2, N-data_exp.N_1, R.B, info, Fsrs_F1(end), Fsrs_del_forward_F1(end), R.kR);
% 
% q_forward  = [q_forward_F1 q_forward_F2];
% qd_forward = [qd_forward_F1 qd_forward_F2];
% lMtilda_forward = [lMtilda_forward_F1 lMtilda_forward_F2];
% Fsrs_del_forward = [Fsrs_del_forward_F1 Fsrs_del_forward_F2];

%% Plot
% Results
h = PlotResults(R,info);
saveas(h,[pathTemp,'/Results/Figures/', info.subj, '_T',num2str(info.trial),info.option,'_1.Results.fig']);

% Params
g = PlotParams(R,info);
saveas(g,[pathTemp,'/Results/Figures/', info.subj, '_T',num2str(info.trial),info.option,'_2.Params.fig']);

% Muscle geometry
% f = PlotMuscleGeometry(R, info);
% saveas(f,[pathTemp,'/Results/Figures/', info.subj, '_T',num2str(info.trial),info.option,'_3.MuscleGeometry.fig']);

% Torques
% p = PlotTorques_SRS(R);
% saveas(f,[pathTemp,'/Results/Figures/', info.subj, '_T',num2str(info.trial),info.option,'_4.Torques.fig']);
