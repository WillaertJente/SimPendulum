%% Simulations of pendulum test using a muscle model - multiple trials together with only kFpe fixed
%  Jente Willaert - 27/04/2021
clear all; close all; clc;
pathmain = pwd;
[pathTemp,~,~] = fileparts(pathmain);
[pathRepo,~,~] = fileparts(pathTemp);

%% Subject and trial information
s.nu   = 'CP6';                                                            % subject number/ name
s.tr   = [4 4];                                                             % subject trials (number of trials)
path   = [pathRepo '\Implicit\Muscle\Experimental data\' s.nu '\'];         % Path to opensim model (scaled)
opt    = '_MT';                                                             % Option used as name to save results
params = ImportParameters(s.nu);                                            % Input parameters (mtot,lc, l, age, m, RG, SE, Nmr, z)

%% Define length of different trials 
% Pre-allocate vectors 
N_all   = zeros(1,length(s.tr)); 
N_1_all = zeros(1,length(s.tr)); 

for j = 1:length(s.tr)
    % Load experimental data for each trial
    [q_exp_r, t_exp, t_span, on_srs] = LoadExpData(s.nu,s.tr,j,0, params, path);        % 1 if you want to plot experimental data
    
    % Discretised time: interpolate experimental data at discr. time using spline
    [q_exp, qdot_exp,N, tvect, dt] = ExpDatAtDiscrTime(t_span,t_exp,q_exp_r);           % dt = 0.005
    
    % Save length of every trial
    N_all(j) = N; 
    
    % Interpolate experimental data 
    xgrid = tvect;
    V     = q_exp;
    W     = qdot_exp;
    q_exp_spline    = casadi.interpolant('LUT','bspline',{xgrid},V);
    q_dot_spline    = casadi.interpolant('LUT','bspline',{xgrid},W);
    
    % Define phases of pendulum (initial state, end of first swing)
    [x0, N_1] = PendulumPhases(q_exp, qdot_exp, N, 1);   
    
    % Save end of first swing for every trial 
    N_1_all(j) = N_1; 
end

%% Define vector for casadi variables with length of all trials together 
import casadi.*;                                                            % Import casadi libraries
opti        = casadi.Opti();                                                % Initialise opti structure
N_tot       = sum(N_all); 
N_tot_Fsrs2 = sum(N_all - N_1_all); 

x_tot            = opti.variable(1,N_tot);                                      % angle (rad)
xd_tot           = opti.variable(1,N_tot);                                      % velocity (rad/s)
lMtilda_ext_tot  = opti.variable(1,N_tot);
lMtilda_flex_tot = opti.variable(1,N_tot);
Fsrs2_tot        = opti.variable(1,N_tot_Fsrs2);                                % SRS force during second phase (exponential decay)
Fsrs_d_tot       = opti.variable(1,N_tot);                                      % Delayed SRS force
a_ext_tot        = opti.variable(1,N_tot);                                      % Activation of extensor

% Controls
lM_projected_ext_tot  = opti.variable(1,N_tot);
lM_projected_flex_tot = opti.variable(1,N_tot);
act_tot               = opti.variable(1,N_tot);                                 % Actuator (non phyisological)
dt1_tot               = opti.variable(length(N_all));

% Slack controls
vMtilda_ext_tot  = opti.variable(1,N_tot);
vMtilda_flex_tot = opti.variable(1,N_tot);

% Parameters
a_ext_0_tot        = opti.variable(length(N_all));                              % Baseline muscle tone extensor
a_flex_tot         = opti.variable(length(N_all));                              % Baseline muscle tone flexor
kFpe_ext_tot       = opti.variable(1);
kFpe_flex_tot      = opti.variable(1);
Rk_tot             = opti.variable(length(N_all));                              % Reflex gain extensor

%% Loop over elke trial (telkens 1 deel van overkoepelende vector selecteren) 
for j = 1:length(s.tr)
    %% Experimental data 
    % Select experimental data for simulation 
    % Load experimental data + SRS on or off
    % Q_exp = BK data - 90°
    [q_exp_r, t_exp, t_span, on_srs] = LoadExpData(s.nu,s.tr,j,0, params, path);      % 1 if you want to plot experimental data
    
    % Discretised time: interpolate experimental data at discr. time using spline
    [q_exp, qdot_exp,N, tvect, dt] = ExpDatAtDiscrTime(t_span,t_exp,q_exp_r);       % dt = 0.005
    
    % Interpolate data 
    xgrid = tvect;
    V     = q_exp;
    W     = qdot_exp;
    q_exp_spline    = casadi.interpolant('LUT','bspline',{xgrid},V);
    q_dot_spline    = casadi.interpolant('LUT','bspline',{xgrid},W);
    
    % Define phases of pendulum (initial state, end of first swing)
    [x0, N_1] = PendulumPhases(q_exp, qdot_exp, N, 1);                          % 1 if you want to plot the phases
    
    %% Information on OpenSim Model
    addpath('MuscleModel');
    
    % Opensim model
    import org.opensim.modeling.*
    model_path = [path,'/',s.nu,'_ScaledModel_ScaledForces_Torso.osim'];                                 % if cp = CPModel_Scaled.osim
    osimModel  = Model(model_path);
    
    % Inertial parameters (tibia)
    % Massa
    bodies         = osimModel.getBodySet();
    if params.z == 18
        tibia      = bodies.get('tibia_l');
        talus      = bodies.get('talus_l');
        calcn      = bodies.get('calcn_l');
        toes       = bodies.get('toes_l');
    else
        tibia      = bodies.get('tibia_r');
        talus      = bodies.get('talus_r');
        calcn      = bodies.get('calcn_r');
        toes       = bodies.get('toes_r');
    end
    params.mass_OStibia = tibia.getMass();
    params.mass_OScalcn = calcn.getMass();
    params.mass_OS = params.mass_OStibia + params.mass_OScalcn;
    
    % COM
    com_tibia = ArrayDouble.createVec3(0);
    tibia.getMassCenter(com_tibia);
    com_calcn = ArrayDouble.createVec3(0);
    calcn.getMassCenter(com_calcn);
    com_foot_in_tibia = abs(com_calcn.get(1)) + params.length_tibia;
    params.lc_OS = (abs(com_tibia.get(1)) * tibia.getMass() + com_foot_in_tibia*calcn.getMass())/params.mass_OS;
    
    % Inertia
    inertia_tibia = Mat33(0);
    tibia.getInertia(inertia_tibia);
    inertia_calcn = Mat33(0);
    calcn.getInertia(inertia_calcn);
    params.I_OS = inertia_tibia.get(0,0) + inertia_calcn.get(0,0) + tibia.getMass() * com_tibia.get(1) ^ 2 + calcn.getMass() * com_foot_in_tibia ^ 2;
    
    % Muscle tendon properties
    [params_Muscle,lOpt,L_TendonSlack,Fiso,PennationAngle]=ReadMuscleParameters(model_path,{'rect_fem_r','bifemlh_r'});
    params_Muscle_ext(2,1) = params_Muscle(2,1);  params_Muscle_flex(2,1) = params_Muscle(2,2);
    params_Muscle_ext(3,1) = params_Muscle(3,1);  params_Muscle_flex(3,1) = params_Muscle(3,2);
    params.MTparams_ext    = params_Muscle(:,1);  params.MTparams_flex    = params_Muscle(:,2); % 5 x 2 matrix: 1 - Fmo, 2 - Lmo, 3 - Lts, 4 - alphao, 5 - vmmax
    
    % Force- velocity characteristics
    load Fvparam.mat
    Fvparam(1)       = 1.475*Fvparam(1); Fvparam(2) = 0.25*Fvparam(2);
    Fvparam(3)       = Fvparam(3) + 0.75; Fvparam(4) = Fvparam(4) - 0.027;
    params.Fvparam   = Fvparam;
    
    % Force length characteristics (active)
    load Faparam.mat
    params.Faparam   = Faparam;
    
    % Force length characteristics (passive)
    e0  = 0.6; kpe = 4; t50  = exp(kpe * (0.2 - 0.10e1) / e0);
    pp1 = (t50 - 0.10e1); t7 = exp(kpe); pp2 = (t7 - 0.10e1);
    params.Fpparam = [pp1;pp2];
    
    %% Initial guess for lM projected en lm tilde vmtilde
    %vMtilde guess
    vMGuess = ones(1,N)*-0.3;
    
    % lMtilde guess
    lMtildeGuess = ones(1,N)*0.9; % 0.9 of 1.1?
    
    % lM projected guess extensor
    lMo_ext      = params.MTparams_ext(2,:)';
    alphao_ext   = params.MTparams_ext(4,:)';
    lMGuess_ext  = lMtildeGuess.*lMo_ext;
    w_ext        = lMo_ext.*sin(alphao_ext);
    lM_projectedGuess_ext = sqrt((lMGuess_ext.^2 - w_ext.^2));
    
    % lM projected guess flexor
    lMo_flex      = params.MTparams_flex(2,:)';
    alphao_flex   = params.MTparams_flex(4,:)';
    lMGuess_flex  = lMtildeGuess.*lMo_flex;
    w_flex        = lMo_flex.*sin(alphao_flex);
    lM_projectedGuess_flex = sqrt((lMGuess_flex.^2 - w_flex.^2));
    
    %dLdMT guess
    dlMdtGuess_ext  = vMGuess*params.MTparams_ext(5)/lMo_ext;
    dlMdtGuess_flex = vMGuess*params.MTparams_flex(5)/lMo_flex;
    
    %% Calculate LMT en Ma
    map_MA         = [path 'MA_FakeMot_T',num2str(s.tr(j))];
    [coeff_LMT_ma_ext, coeff_LMT_ma_flex] = DefineLMTCoefficients_2muscles_MT(map_MA, s.nu);
    
    %% Calculate offset (difference between IK en BK)
    % BK + offset = IK
    IKtrial  = [path,'IK_Trial0',num2str(s.tr(j)),'.mot'];
    [offset] = CalculateOffset(IKtrial, params);
    
    %% Define states, controls, bounds and initial guess for every trial seperate     
    
    % States
    % xtrial = x(1:N1)
    if j == 1
        x1             = x_tot(1,1:N_all(j));                             % angle (rad)
        xd1            = xd_tot(1,1:N_all(j));                             % velocity (rad/s)
        lMtilda_ext1   = lMtilda_ext_tot(1,1:N_all(j));
        lMtilda_flex1  = lMtilda_flex_tot(1,1:N_all(j));
        Fsrs21         = Fsrs2_tot(1,1:N_all(j)-N_1_all(j));      % SRS force during second phase (exponential decay)
        Fsrs_d1        = Fsrs_d_tot(1,1:N_all(j));                             % Delayed SRS force
        a_ext1         = a_ext_tot(1,1:N_all(j));                             % Activation of extensor
        
        % Controls
        lM_projected_ext1  = lM_projected_ext_tot(1,1:N_all(j));
        lM_projected_flex1 = lM_projected_flex_tot(1,1:N_all(j));
        act1               = act_tot(1,1:N_all(j));     % Actuator (non phyisological)
        dt11               = dt1_tot(1);
        
        % Slack controls
        vMtilda_ext1  = vMtilda_ext_tot(1,1:N_all(j));
        vMtilda_flex1 = vMtilda_flex_tot(1,1:N_all(j));
        
        % Parameters
        a_ext_01        = a_ext_0_tot(j);                                    % Baseline muscle tone extensor
        a_flex1         = a_flex_tot(j);                                     % Baseline muscle tone flexor
        kFpe_ext1       = kFpe_ext_tot;
        kFpe_flex1      = kFpe_flex_tot;
        Rk1             = Rk_tot(j);                                         % Reflex gain extensor
        
        % Bounds
    opti.subject_to(-4*pi < x1   < 4*pi);        % Friedl - It seems that these bounds are assuming angles are in degrees, not the case. Smaller/larger bounds => no convergence.
    opti.subject_to(-300  < xd1  < 300);
    opti.subject_to(0.0   < a_ext_01 < 0.5);
    opti.subject_to(0.0   < a_flex1  < 0.5);
    opti.subject_to(0.0   < a_ext1   < 0.5);
    opti.subject_to(1e-4  < lM_projected_ext1);  % Only positive lM's
    opti.subject_to(1e-4  < lM_projected_flex1);
    opti.subject_to(-10   < vMtilda_ext1 < 10);  % Friedl - changed bounds
    opti.subject_to(-10   < vMtilda_flex1 < 10); % Friedl - changed bounds
    opti.subject_to(0.2   < lMtilda_ext1 < 1.8);
    opti.subject_to(0.2   < lMtilda_flex1 < 1.8);
    opti.subject_to(0.0   < Fsrs21   < 2);
    opti.subject_to(-1    < Fsrs_d1  < 2);
    opti.subject_to(0.055 < kFpe_ext1    < 0.15);
    opti.subject_to(0.05  < kFpe_flex1   < 0.15);
    opti.subject_to(0.001< dt11     < 0.01);    % 0.05
    opti.subject_to(1e-5    < Rk1      < 4);
%     opti.subject_to(-0.005< act < 0.005);
    
    % Bounds on initial states
    opti.subject_to(x1(1)     == x0(1));
    opti.subject_to(xd1(1)    == 0);
    opti.subject_to(Fsrs_d1(1) ==0);
    
    % Initial guess
    opti.set_initial(x1, q_exp);
    opti.set_initial(xd1,qdot_exp);
    opti.set_initial(a_ext_01, 0.001);
    opti.set_initial(a_flex1,  0.01);
    opti.set_initial(a_ext1,   0.001*ones(1,N));
    opti.set_initial(kFpe_ext1,0.1);
    opti.set_initial(kFpe_flex1,0.1);
    opti.set_initial(lM_projected_ext1, lM_projectedGuess_ext);
    opti.set_initial(lM_projected_flex1, lM_projectedGuess_flex);
    opti.set_initial(lMtilda_ext1, lMtildeGuess);
    opti.set_initial(lMtilda_flex1, lMtildeGuess);
    opti.set_initial(vMtilda_ext1, vMGuess);
    opti.set_initial(vMtilda_flex1, vMGuess);
    opti.set_initial(act1,0);
    opti.set_initial(dt11,0.005);
    opti.set_initial(Rk1, 0.01);
    
    % Defining problem (muscle model)
    % Calculate shift
    kT = 35;
    [shift]    = getshift(kT);
    
    offset   =  mean(offset)*pi/180; % Friedl - moved this, you do not need to compute this every time you call muscle dynamics.
    
    % Calculate FT en ma
    %[FT_ext,FT_flex, ma_ext, ma_flex, dlMdt_ext, dlMdt_flex, err_ext, err_flex, lM_ext, lM_flex, lT_ext, lT_flex, Fce_ext, Fce_flex, Fpe_ext, Fpe_flex, FM_ext, FM_flex, Fsrs, Fsrs_dot, FMltilda_ext, FMltilda_flex] = CalculateTendonForceAndMomentArm_v3_2muscles(x, params, lMtilda_ext, lMtilda_flex, a_ext,a_flex, shift, vMtilda_ext,vMtilda_flex, lM_projected_ext, lM_projected_flex, coeff_LMT_ma_ext, coeff_LMT_ma_flex, offset, kFpe_ext, kFpe_flex, N_1,Fsrs, N);
    [FT_ext1,FT_flex1, ma_ext1, ma_flex1, dlMdt_ext1, dlMdt_flex1, err_ext1, err_flex1, ...
        lM_ext1, lM_flex1, lT_ext1, lT_flex1, Fce_ext1, Fce_flex1, Fpe_ext1, Fpe_flex1, ...
        FM_ext1, FM_flex1, Fsrs11, FMltilda_ext1, FMltilda_flex1, FMvtilda_ext1] = ...
        CalculateTendonForceAndMomentArm_v7_Fmv(x1, params, lMtilda_ext1, ...
        lMtilda_flex1, a_ext_01, a_flex1, shift, vMtilda_ext1, vMtilda_flex1, ...
        lM_projected_ext1, lM_projected_flex1, coeff_LMT_ma_ext, coeff_LMT_ma_flex, ...
        offset, kFpe_ext1, kFpe_flex1, N_1, N, Fsrs_d1, Fsrs21, a_ext1);
    
    opti.subject_to(-0.0001 < a_ext1 - a_ext_01 - Rk1*Fsrs_d1 <0.0001); % Does not converge with smaller bounds.
    %     opti.subject_to(a_ext == a_ext_0 + Rk*Fsrs_d); % Friedl - does not
    %     converge. Still wondering why but not surprising given our previous
    %     results. Also, when eliminating a_ext from the problem, it did not
    %     converge.
    
    % Variable time of phases
    tF11 = (N_1)*dt11;        %N_1-1
    dt2 = 0.005;
    
    % Constraints on phases
    opti.subject_to(tF11 > 0.2 );
    opti.subject_to(xd1(1:N_1) < 0);
    opti.subject_to(1e-4 < xd1(N_1+1) < 1); % positive joint velocity after end of first swing
    
    % Constraints srs - exponential decay
    opti.subject_to(Fsrs21 == [Fsrs11(N_1) Fsrs21(1:N-N_1-1)] - [Fsrs11(N_1) Fsrs21(1:N-N_1-1)]/ 0.05* dt2); % Friedl
    
    % Reflex activity
    tau_d = 0.05;
    Fsrs_ddt1          = [(Fsrs11-Fsrs_d1(1:N_1))  (Fsrs21-Fsrs_d1(N_1+1:N))]/tau_d;
    
    % Dynamics
    xdd1 = 1/params.I_OS * ((-params.mass_OS*params.g*params.lc_OS*cos(x1))+ FT_ext1.*ma_ext1 + FT_flex1.*ma_flex1 + act1 - 0.0117*xd1); %  + Tdamp + FT*ma);
    
    % backward euler
    % opti.subject_to(xd(1:N-1)*dt +x(1:N-1) == x(2:N));
    opti.subject_to(xd1(1:N_1-1)*dt11 +x1(1:N_1-1) == x1(2:N_1));
    opti.subject_to(xd1(N_1:N-1)*dt2 +x1(N_1:N-1) == x1(N_1+1:N));
    % opti.subject_to(xdd(1:N-1)*dt +xd(1:N-1) == xd(2:N));
    opti.subject_to(xdd1(1:N_1-1)*dt11 +xd1(1:N_1-1) == xd1(2:N_1));
    opti.subject_to(xdd1(N_1:N-1)*dt2 +xd1(N_1:N-1) == xd1(N_1+1:N));
    % opti.subject_to(dlMdt_ext(1:N-1)*dt + lMtilda_ext(1:N-1) == lMtilda_ext(2:N)); % VmTilde met factor 10 (MRS)
    opti.subject_to(dlMdt_ext1(1:N_1-1)*dt11 + lMtilda_ext1(1:N_1-1) == lMtilda_ext1(2:N_1));
    opti.subject_to(dlMdt_ext1(N_1:N-1)*dt2 + lMtilda_ext1(N_1:N-1) == lMtilda_ext1(N_1+1:N));
    % opti.subject_to(dlMdt_flex(1:N-1)*dt + lMtilda_flex(1:N-1) == lMtilda_flex(2:N)); % VmTilde met factor 10 (MRS)
    opti.subject_to(dlMdt_flex1(1:N_1-1)*dt11 + lMtilda_flex1(1:N_1-1) == lMtilda_flex1(2:N_1));
    opti.subject_to(dlMdt_flex1(N_1:N-1)*dt2 + lMtilda_flex1(N_1:N-1) == lMtilda_flex1(N_1+1:N));
    opti.subject_to(Fsrs_ddt1(1:N_1-1)*dt11 + Fsrs_d1(1:N_1-1) == Fsrs_d1(2:N_1));
    opti.subject_to(Fsrs_ddt1(N_1:N-1)*dt2 + Fsrs_d1(N_1:N-1) == Fsrs_d1(N_1+1:N));
    opti.subject_to(lM_ext1.^2 - w_ext.^2 == lM_projected_ext1.^2);
    opti.subject_to(lM_flex1.^2 - w_flex.^2 == lM_projected_flex1.^2);
    opti.subject_to(err_ext1 == 0);
    opti.subject_to(err_flex1 == 0);
    
    % Objective function
    tvect_spline = MX(1,N);
    for i = 1:N_1
        t = tvect(1)+dt11*i;
        tvect_spline(i) = t;
    end
    
    for i = N_1+1:N
        t = tvect(1)+dt11*N_1+dt2*i;
        tvect_spline(i) = t;
    end
    
    error        = x1 - q_exp;
    error_dot    = xd1 - qdot_exp;
    
    J1            = sumsqr(error)  + sumsqr(error_dot)  + 100*sumsqr(act1)+ 0.001 * (sumsqr(vMtilda_ext1))  %+ 100*sumsqr(act);%+ ...
    % 0.001 * (sumsqr(vMtilda_ext)); %+ 100*sumsqr(error_ra);
    opti.minimize(J1);
        
    else  
        x2             = x_tot(1,N_all(j-1)+1:N_all(j-1)+N_all(j));
        xd2            = xd_tot(1,N_all(j-1)+1:N_all(j-1)+N_all(j));
        lMtilda_ext2   = lMtilda_ext_tot(1,N_all(j-1)+1:N_all(j-1)+N_all(j));
        lMtilda_flex2  = lMtilda_flex_tot(1,N_all(j-1)+1:N_all(j-1)+N_all(j));
        Fsrs22         = Fsrs2_tot(N_all(j-1)-N_1_all(j-1)+1:N_all(j-1)-N_1_all(j-1)+N_all(j)-N_1_all(j));
        Fsrs_d2        = Fsrs_d_tot(1,N_all(j-1)+1:N_all(j-1)+N_all(j)); 
        a_ext2         = a_ext_tot(1,N_all(j-1)+1:N_all(j-1)+N_all(j));    
        
        % Controls
        lM_projected_ext2  =  lM_projected_ext_tot(1,N_all(j-1)+1:N_all(j-1)+N_all(j));
        lM_projected_flex2 = lM_projected_flex_tot(1,N_all(j-1)+1:N_all(j-1)+N_all(j));
        act2               = act_tot(1,N_all(j-1)+1:N_all(j-1)+N_all(j));     % Actuator (non phyisological)
        dt12               = dt1_tot(j);
        
        % Slack controls
        vMtilda_ext2  = vMtilda_ext_tot(1,N_all(j-1)+1:N_all(j-1)+N_all(j));
        vMtilda_flex2 = vMtilda_flex_tot(1,N_all(j-1)+1:N_all(j-1)+N_all(j));
        
        % Parameters
        a_ext_02        = a_ext_0_tot(j);                                    % Baseline muscle tone extensor
        a_flex2         = a_flex_tot(j);                                     % Baseline muscle tone flexor
        kFpe_ext2       = kFpe_ext_tot;
        kFpe_flex2      = kFpe_flex_tot;
        Rk2             = Rk_tot(j);   

    % Bounds
    opti.subject_to(-4*pi < x2   < 4*pi);        % Friedl - It seems that these bounds are assuming angles are in degrees, not the case. Smaller/larger bounds => no convergence.
    opti.subject_to(-300  < xd2  < 300);
    opti.subject_to(0.0   < a_ext_02 < 0.5);
    opti.subject_to(0.0   < a_flex2  < 0.5);
    opti.subject_to(0.0   < a_ext2   < 0.5);
    opti.subject_to(1e-4  < lM_projected_ext2);  % Only positive lM's
    opti.subject_to(1e-4  < lM_projected_flex2);
    opti.subject_to(-10   < vMtilda_ext2 < 10);  % Friedl - changed bounds
    opti.subject_to(-10   < vMtilda_flex2 < 10); % Friedl - changed bounds
    opti.subject_to(0.2   < lMtilda_ext2 < 1.8);
    opti.subject_to(0.2   < lMtilda_flex2 < 1.8);
    opti.subject_to(0.0   < Fsrs22   < 2);
    opti.subject_to(-1    < Fsrs_d2  < 2);
    opti.subject_to(0.055 < kFpe_ext2   < 0.15);
    opti.subject_to(0.05  < kFpe_flex2   < 0.15);
    opti.subject_to(0.001< dt12     < 0.01);    % 0.05
    opti.subject_to(1e-5    < Rk2      < 4);
%     opti.subject_to(-0.005< act < 0.005);
    
    % Bounds on initial states
    opti.subject_to(x2(1)     == x0(1));
    opti.subject_to(xd2(1)    == 0);
    opti.subject_to(Fsrs_d2(1) ==0);
    
    % Initial guess
    opti.set_initial(x2, q_exp);
    opti.set_initial(xd2,qdot_exp);
    opti.set_initial(a_ext_02, 0.001);
    opti.set_initial(a_flex2,  0.01);
    opti.set_initial(a_ext2,   0.001*ones(1,N));
    opti.set_initial(kFpe_ext2,0.1);
    opti.set_initial(kFpe_flex2,0.1);
    opti.set_initial(lM_projected_ext2, lM_projectedGuess_ext);
    opti.set_initial(lM_projected_flex2, lM_projectedGuess_flex);
    opti.set_initial(lMtilda_ext2, lMtildeGuess);
    opti.set_initial(lMtilda_flex2, lMtildeGuess);
    opti.set_initial(vMtilda_ext2, vMGuess);
    opti.set_initial(vMtilda_flex2, vMGuess);
    opti.set_initial(act2,0);
    opti.set_initial(dt12,0.005);
    opti.set_initial(Rk2, 0.01);
    
    % Defining problem (muscle model)
    % Calculate shift
    kT = 35;
    [shift]    = getshift(kT);
    
    offset   =  mean(offset)*pi/180; % Friedl - moved this, you do not need to compute this every time you call muscle dynamics.
    
    % Calculate FT en ma
    %[FT_ext,FT_flex, ma_ext, ma_flex, dlMdt_ext, dlMdt_flex, err_ext, err_flex, lM_ext, lM_flex, lT_ext, lT_flex, Fce_ext, Fce_flex, Fpe_ext, Fpe_flex, FM_ext, FM_flex, Fsrs, Fsrs_dot, FMltilda_ext, FMltilda_flex] = CalculateTendonForceAndMomentArm_v3_2muscles(x, params, lMtilda_ext, lMtilda_flex, a_ext,a_flex, shift, vMtilda_ext,vMtilda_flex, lM_projected_ext, lM_projected_flex, coeff_LMT_ma_ext, coeff_LMT_ma_flex, offset, kFpe_ext, kFpe_flex, N_1,Fsrs, N);
    [FT_ext2,FT_flex2, ma_ext2, ma_flex2, dlMdt_ext2, dlMdt_flex2, err_ext2, err_flex2, ...
        lM_ext2, lM_flex2, lT_ext2, lT_flex2, Fce_ext2, Fce_flex2, Fpe_ext2, Fpe_flex2, ...
        FM_ext2, FM_flex2, Fsrs12, FMltilda_ext2, FMltilda_flex2, FMvtilda_ext2] = ...
        CalculateTendonForceAndMomentArm_v7_Fmv(x2, params, lMtilda_ext2, ...
        lMtilda_flex2, a_ext_02, a_flex2, shift, vMtilda_ext2, vMtilda_flex2, ...
        lM_projected_ext2, lM_projected_flex2, coeff_LMT_ma_ext, coeff_LMT_ma_flex, ...
        offset, kFpe_ext2, kFpe_flex2, N_1, N, Fsrs_d2, Fsrs22, a_ext2);
    
    opti.subject_to(-0.0001 < a_ext2 - a_ext_02 - Rk2*Fsrs_d2 <0.0001); % Does not converge with smaller bounds.
    %     opti.subject_to(a_ext == a_ext_0 + Rk*Fsrs_d); % Friedl - does not
    %     converge. Still wondering why but not surprising given our previous
    %     results. Also, when eliminating a_ext from the problem, it did not
    %     converge.
    
    % Variable time of phases
    tF1 = (N_1)*dt12;        %N_1-1
    dt2 = 0.005;
    
    % Constraints on phases
    opti.subject_to(tF1 > 0.2 );
    opti.subject_to(xd2(1:N_1) < 0);
    opti.subject_to(1e-4 < xd2(N_1+1) < 1); % positive joint velocity after end of first swing
    
    % Constraints srs - exponential decay
    opti.subject_to(Fsrs22 == [Fsrs12(N_1) Fsrs22(1:N-N_1-1)] - [Fsrs12(N_1) Fsrs22(1:N-N_1-1)]/ 0.05* dt2); % Friedl
    
    % Reflex activity
    tau_d = 0.05;
    Fsrs_ddt2          = [(Fsrs12-Fsrs_d2(1:N_1))  (Fsrs22-Fsrs_d2(N_1+1:N))]/tau_d;
    
    % Dynamics
    xdd2 = 1/params.I_OS * ((-params.mass_OS*params.g*params.lc_OS*cos(x2))+ FT_ext2.*ma_ext2 + FT_flex2.*ma_flex2 + act2 - 0.0117*xd2); %  + Tdamp + FT*ma);
    
    % backward euler
    % opti.subject_to(xd(1:N-1)*dt +x(1:N-1) == x(2:N));
    opti.subject_to(xd2(1:N_1-1)*dt12 +x2(1:N_1-1) == x2(2:N_1));
    opti.subject_to(xd2(N_1:N-1)*dt2 +x2(N_1:N-1) == x2(N_1+1:N));
    % opti.subject_to(xdd(1:N-1)*dt +xd(1:N-1) == xd(2:N));
    opti.subject_to(xdd2(1:N_1-1)*dt12 +xd2(1:N_1-1) == xd2(2:N_1));
    opti.subject_to(xdd2(N_1:N-1)*dt2 +xd2(N_1:N-1) == xd2(N_1+1:N));
    % opti.subject_to(dlMdt_ext(1:N-1)*dt + lMtilda_ext(1:N-1) == lMtilda_ext(2:N)); % VmTilde met factor 10 (MRS)
    opti.subject_to(dlMdt_ext2(1:N_1-1)*dt12 + lMtilda_ext2(1:N_1-1) == lMtilda_ext2(2:N_1));
    opti.subject_to(dlMdt_ext2(N_1:N-1)*dt2 + lMtilda_ext2(N_1:N-1) == lMtilda_ext2(N_1+1:N));
    % opti.subject_to(dlMdt_flex(1:N-1)*dt + lMtilda_flex(1:N-1) == lMtilda_flex(2:N)); % VmTilde met factor 10 (MRS)
    opti.subject_to(dlMdt_flex2(1:N_1-1)*dt12 + lMtilda_flex2(1:N_1-1) == lMtilda_flex2(2:N_1));
    opti.subject_to(dlMdt_flex2(N_1:N-1)*dt2 + lMtilda_flex2(N_1:N-1) == lMtilda_flex2(N_1+1:N));
    opti.subject_to(Fsrs_ddt2(1:N_1-1)*dt12 + Fsrs_d2(1:N_1-1) == Fsrs_d2(2:N_1));
    opti.subject_to(Fsrs_ddt2(N_1:N-1)*dt2 + Fsrs_d2(N_1:N-1) == Fsrs_d2(N_1+1:N));
    opti.subject_to(lM_ext2.^2 - w_ext.^2 == lM_projected_ext2.^2);
    opti.subject_to(lM_flex2.^2 - w_flex.^2 == lM_projected_flex2.^2);
    opti.subject_to(err_ext2 == 0);
    opti.subject_to(err_flex2 == 0);
    
    % Objective function
    tvect_spline = MX(1,N);
    for i = 1:N_1
        t = tvect(1)+dt12*i;
        tvect_spline(i) = t;
    end
    
    for i = N_1+1:N
        t = tvect(1)+dt12*N_1+dt2*i;
        tvect_spline(i) = t;
    end
    
    error        = x2 - q_exp;
    error_dot    = xd2 - qdot_exp;
    
    J2            = sumsqr(error)  + sumsqr(error_dot)  + 100*sumsqr(act2)+ 0.001 * (sumsqr(vMtilda_ext2)) %+ 100*sumsqr(act);%+ ...
    % 0.001 * (sumsqr(vMtilda_ext)); %+ 100*sumsqr(error_ra);
    opti.minimize(J2);
end
end
    % options for IPOPT
    options.ipopt.tol = 1*10^(-6);          % moet normaal 10^-6 zijn
    options.ipopt.linear_solver = 'mumps';
    %     options.ipopt.nlp_scaling_method = 'none';
    %     options.ipopt.linear_solver = 'ma57';
    %options.ipopt.hessian_approximation = 'limited-memory'; % enkel bij moeilijke problemen
    
    % Solve the OCP
    opti.solver('ipopt',options);
    sol = opti.solve();
    
    sol_x = sol.value(x_tot);                       sol_dx = sol.value(xd_tot);
    sol_a_ext = sol.value(a_ext_tot);               sol_a_flex = sol.value(a_flex_tot);
    sol_lMtilda_ext = sol.value(lMtilda_ext_tot);   sol_lMtilda_flex = sol.value(lMtilda_flex_tot);
    sol_aext0 = sol.value(a_ext_0_tot);
    sol_act = sol.value(act_tot);
    sol_FT_ext  = sol.value(FT_ext_tot);            sol_FT_flex  = sol.value(FT_flex_tot);
    sol_ma_ext = sol.value(ma_ext_tot);             sol_ma_flex = sol.value(ma_flex_tot);
    sol_Fpe_ext = sol.value(Fpe_ext_tot);           sol_Fpe_flex = sol.value(Fpe_flex_tot);
    %sol_kFpe    = sol.value(kFpe);
    sol_kFpe_ext = sol.value(kFpe_ext_tot);
    sol_kFpe_flex = sol.value(kFpe_flex_tot);
    %sol_kFpe_ext = sol.value(kFpe_ext);         sol_kFpe_flex = sol.value(kFpe_flex);
%     sol_J   = sol.value(J);
%     sol_Fsrs1 = sol.value(Fsrs1_tot);
    sol_Fsrs2 = sol.value(Fsrs2_tot);
    %sol_Fsrs = [sol_Fsrs1 sol_Fsrs2];
    sol_Fsrs_d = sol.value(Fsrs_d_tot);
    sol_dt1  = sol.value(dt1_tot);
    sol_Rk   = sol.value(Rk_tot);
    %sol_dlmdt_ext = sol.value(dlMdt_ext_tot);       sol_dlmdt_flex = sol.value(dlMdt_flex);
    %sol_lM_ext = sol.value(lM_ext_tot);             sol_lM_flex    = sol.value(lM_flex);
    %sol_FMltilda_ext = sol.value(FMltilda_ext_tot); sol_FMltilda_flex = sol.value(FMltilda_flex);
    %Total = sol.value(FT_ext).*sol.value(ma_ext) + sol.value(FT_flex).*sol.value(ma_flex);
    save(['C:\Users\u0125183\Box\PhD 1\Simulations Pendulum Test\Results/Result_',char(s.nu),char(opt),'.mat'],'sol_x','sol_dx','sol_a_ext','sol_a_flex','sol_aext0', 'sol_a_flex','sol_a_ext', 'sol_lMtilda_ext','sol_lMtilda_flex', 'sol_FT_ext', 'sol_FT_flex', 'sol_ma_ext', 'sol_ma_flex', 'sol_Fpe_ext', 'sol_Fpe_flex', 'sol_kFpe_ext','sol_kFpe_flex', 'sol_J', 'sol_Fsrs','q_exp','sol_dlmdt_ext','sol_dlmdt_flex','sol_lM_ext','sol_lM_flex','sol_FMltilda_ext','sol_FMltilda_flex','tvect','sol_Fsrs_d')
    
    figure(j*10)
    plot(tvect,q_exp,'k','LineWidth',1.5)
    hold on
    plot(tvect,sol_x,'LineWidth',1.5)
    
    figure(j*100)
    subplot(611)
    plot(tvect,q_exp*180/pi,'k','LineWidth',1.5)
    hold on
    plot(tvect,sol_x*180/pi,'LineWidth',1.5)
    hold on
    ylabel('({\circ})'); box off; legend('Exp','Sim');
    subplot(612)
    plot(tvect,sol.value(FT_ext),'LineWidth',1.5);
    hold on; plot(tvect,sol.value(FT_flex),'LineWidth',1.5);
    hold on; box off; ylabel('FT');
    %     subplot(613)
    %     plot(tvect,sol.value(ma_ext),'LineWidth',1.5);
    %     hold on; plot(tvect,sol.value(ma_flex),'LineWidth',1.5);
    %     hold on; box off; ylabel('ma');
    subplot(613)
    plot(tvect,sol.value(Fpe_ext),'LineWidth',1.5);
    hold on;  plot(tvect,sol.value(Fpe_flex),'LineWidth',1.5);
    hold on; box off; ylabel('Fpe');
    subplot(614)
    plot(tvect,sol_Fsrs,'LineWidth',1.5); hold on;
    plot(tvect,sol_Fsrs_d,'LineWidth',1); box off; ylabel('Fsrs');
    subplot(615)
    plot(tvect,sol_aext0 + sol.value(Rk)*sol_Fsrs_d,'LineWidth',1.5);
    hold on; box off; ylabel('a0+Rk*Fsrs_d');
    subplot(616)
    plot(tvect,sol.value(act),'LineWidth',1.5);
    hold on;
    plot(tvect, Total,'LineWidth',1.5);
    hold on; box off; ylabel('Nm');  legend({'actuator','Total torque muscles'})
    
    figure()
    plot(sol.value(act),'LineWidth',1.5);
    disp(['Reflex gain =', num2str(sol.value(Rk))])
    
    figure()
    bar(1,sol_kFpe_ext); hold on;
    bar(2,sol_kFpe_flex); hold on;
    legend('Ext','Flex');

