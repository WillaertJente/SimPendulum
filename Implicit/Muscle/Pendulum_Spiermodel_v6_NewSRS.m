%% Simulations of pendulum test using a muscle model
%  Jente Willaert - 21/10/2020
% clear all; close all; clc;

%% Input
s.nu = 'CP7';                                                                               % subject number/ name
s.tr = [2];                                                                                % subject trials (number of trials)
pathmain = pwd;
[pathTemp,~,~] = fileparts(pathmain);
[pathRepo,~,~] = fileparts(pathTemp);
path = [pathRepo '\Implicit\Muscle\Experimental data\' s.nu '\'];                       % Path to opensim model (scaled)
opt  = '_NewSRS_N1_30';   % Option used as name to save results

params = ImportParameters(s.nu);    % Input parameters (mtot,lc, l, age, m, RG, SE, Nmr, z)

for j = 1:length(s.tr)
    
    %% Prepare experimental data for simulation
    % Load experimental data + SRS on or off
    % Q_exp = BK data - 90°
    [q_exp_r, t_exp, t_span, on_srs] = LoadExpData(s.nu,s.tr,j,0, params, path);      % 1 if you want to plot experimental data
    
    % Discretised time: interpolate experimental data at discr. time using spline
    [q_exp, qdot_exp,N, tvect, dt] = ExpDatAtDiscrTime(t_span,t_exp,q_exp_r);       % dt = 0.005
    
%       q_exp = q_exp(1:end-300);
%       qdot_exp = qdot_exp(1:end-300);
%       N = N-300;
%       tvect = tvect(1:end-300);
   
    %% Formulate OCP
    import casadi.*;        % Import casadi libraries
    opti  = casadi.Opti();   % Initialise opti structure
    xgrid = tvect;
    V     = q_exp;
    W     = qdot_exp;
    q_exp_spline    = casadi.interpolant('LUT','bspline',{xgrid},V);
    q_dot_spline    = casadi.interpolant('LUT','bspline',{xgrid},W);
    
    % Define phases of pendulum (initial state, end of first swing)
    [x0, N_1] = PendulumPhases(q_exp, qdot_exp, N, 1);                          % 1 if you want to plot the phases
    N_1 = 30; 
    %% Add OpenSim Model
    addpath('MuscleModel');
    
    % Opensim model
    import org.opensim.modeling.*
    model_path = [path,'/',s.nu,'_ScaledModel_ScaledForces.osim'];                                 % if cp = CPModel_Scaled.osim
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
    [params_Muscle,lOpt,L_TendonSlack,Fiso,PennationAngle]=ReadMuscleParameters(model_path,{'rect_fem_l','bifemlh_l'});
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
    [coeff_LMT_ma_ext, coeff_LMT_ma_flex] = DefineLMTCoefficients_2muscles(map_MA, s.nu);
    
    %% Calculate offset (difference between IK en BK)
    % BK + offset = IK
    IKtrial  = [path,'IK_Trial0',num2str(s.tr(j)),'.mot'];
    [offset] = CalculateOffset(IKtrial, params);
    
    %% Define states, controls, bounds and initial guess
    
    % States
    x            = opti.variable(1,N);          % angle (rad)
    xd           = opti.variable(1,N);          % velocity (rad/s)
    lMtilda_ext  = opti.variable(1,N);          
    lMtilda_flex = opti.variable(1,N);
    Fsrs2        = opti.variable(1,N-N_1);      % SRS force during second phase (exponential decay)
    Fsrs_d       = opti.variable(1,N);          % Delayed SRS force
    a_ext        = opti.variable(1,N);          % Activation of extensor 
    
    % Controls
    lM_projected_ext  = opti.variable(1,N);
    lM_projected_flex = opti.variable(1,N);
    act               = opti.variable(1,N);     % Actuator (non phyisological)
    dt1               = opti.variable(1);
    
    % Slack controls
    vMtilda_ext  = opti.variable(1,N);
    vMtilda_flex = opti.variable(1,N);
    
    % Parameters
    a_ext_0        = opti.variable(1);          % Baseline muscle tone extensor
    a_flex         = opti.variable(1);          % Baseline muscle tone flexor
    kFpe_ext       = opti.variable(1); 
    kFpe_flex      = opti.variable(1); 
    Rk             = opti.variable(1);          % Reflex gain extensor 
    
    % Bounds
    opti.subject_to(-4*pi < x   < 4*pi);        % Friedl - It seems that these bounds are assuming angles are in degrees, not the case. Smaller/larger bounds => no convergence.
    opti.subject_to(-300  < xd  < 300);
    opti.subject_to(0.0   < a_ext_0 < 0.5);
    opti.subject_to(0.0   < a_flex  < 0.5);
    opti.subject_to(0.0   < a_ext   < 0.5);
    opti.subject_to(1e-4  < lM_projected_ext);  % Only positive lM's
    opti.subject_to(1e-4  < lM_projected_flex);
    opti.subject_to(-10   < vMtilda_ext  < 10);  % Friedl - changed bounds
    opti.subject_to(-10   < vMtilda_flex < 10); % Friedl - changed bounds
    opti.subject_to(0.2   < lMtilda_ext  < 1.8);
    opti.subject_to(0.2   < lMtilda_flex < 1.8);
    opti.subject_to(0.0   < Fsrs2   < 2);
    opti.subject_to(-1    < Fsrs_d  < 2);
    opti.subject_to(0.05  < kFpe_ext    < 0.15);
    opti.subject_to(0.05  < kFpe_flex   < 0.15);
    opti.subject_to(0.001 < dt1     < 0.01);    % 0.05
    opti.subject_to(1e-5  < Rk      < 4);
    opti.subject_to(-0.01 < act < 0.01); 
    
    % Bounds on initial states
    opti.subject_to(x(1)      == x0(1));
    opti.subject_to(xd(1)     == 0);
    opti.subject_to(Fsrs_d(1) == 0);
    
    % Initial guess
    opti.set_initial(x, q_exp);
    opti.set_initial(xd,qdot_exp);
    opti.set_initial(a_ext_0, 0.001);
    opti.set_initial(a_flex,  0.01);
    opti.set_initial(a_ext,   0.001*ones(1,N));
    opti.set_initial(kFpe_ext,0.1);
    opti.set_initial(kFpe_flex,0.1);
    opti.set_initial(lM_projected_ext,  lM_projectedGuess_ext);
    opti.set_initial(lM_projected_flex, lM_projectedGuess_flex);
    opti.set_initial(lMtilda_ext,  lMtildeGuess);        
    opti.set_initial(lMtilda_flex, lMtildeGuess);       
    opti.set_initial(vMtilda_ext,  vMGuess);
    opti.set_initial(vMtilda_flex, vMGuess);
    opti.set_initial(act,0);
    opti.set_initial(dt1,0.0005);
    opti.set_initial(Rk, 0.01);
    
    % Defining problem (muscle model)
    % Calculate shift
    kT = 35;
    [shift]    = getshift(kT);
    
    offset   =  mean(offset)*pi/180; % Friedl - moved this, you do not need to compute this every time you call muscle dynamics.
    
    % Calculate FT en ma
    %[FT_ext,FT_flex, ma_ext, ma_flex, dlMdt_ext, dlMdt_flex, err_ext, err_flex, lM_ext, lM_flex, lT_ext, lT_flex, Fce_ext, Fce_flex, Fpe_ext, Fpe_flex, FM_ext, FM_flex, Fsrs, Fsrs_dot, FMltilda_ext, FMltilda_flex] = CalculateTendonForceAndMomentArm_v3_2muscles(x, params, lMtilda_ext, lMtilda_flex, a_ext,a_flex, shift, vMtilda_ext,vMtilda_flex, lM_projected_ext, lM_projected_flex, coeff_LMT_ma_ext, coeff_LMT_ma_flex, offset, kFpe_ext, kFpe_flex, N_1,Fsrs, N);
    [FT_ext,FT_flex, ma_ext, ma_flex, dlMdt_ext, dlMdt_flex, err_ext, err_flex, ...
        lM_ext, lM_flex, lT_ext, lT_flex, Fce_ext, Fce_flex, Fpe_ext, Fpe_flex, ...
        FM_ext, FM_flex, Fsrs1, FMltilda_ext, FMltilda_flex] = ...
        CalculateTendonForceAndMomentArm_v4_NewSRS(x, params, lMtilda_ext, ...
        lMtilda_flex, a_ext_0, a_flex, shift, vMtilda_ext, vMtilda_flex, ...
        lM_projected_ext, lM_projected_flex, coeff_LMT_ma_ext, coeff_LMT_ma_flex, ...
        offset, kFpe_ext, kFpe_flex, N_1, N, Fsrs_d, Fsrs2, a_ext);
    
    opti.subject_to(-0.0001 < a_ext - a_ext_0 - Rk*Fsrs_d <0.0001); % Does not converge with smaller bounds.
    %     opti.subject_to(a_ext == a_ext_0 + Rk*Fsrs_d); % Friedl - does not
    %     converge. Still wondering why but not surprising given our previous
    %     results. Also, when eliminating a_ext from the problem, it did not
    %     converge.
    
    % Variable time of phases
    tF1 = (N_1)*dt1;        %N_1-1
    dt2 = 0.005;
    
    % Constraints on phases
    opti.subject_to(tF1 > 0.02 );
    opti.subject_to(lMtilda_ext(1:N_1) - lMtilda_ext(1) <= 5.7*10^-3); 
    opti.subject_to(lMtilda_ext(N_1+1) - lMtilda_ext(1) > 5.7*10^-3);
    
    opti.subject_to(xd(1:N_1) < 0);
%     opti.subject_to(1e-4 < xd(N_1+1) < 1); % positive joint velocity after end of first swing
    
    % Constraints srs - exponential decay
    opti.subject_to(Fsrs2 == [Fsrs1(N_1) Fsrs2(1:N-N_1-1)] - [Fsrs1(N_1) Fsrs2(1:N-N_1-1)]/ 0.05* dt2); % Friedl
    
    % Reflex activity
    tau_d = 0.05;
    Fsrs_ddt          = [(Fsrs1-Fsrs_d(1:N_1))  (Fsrs2-Fsrs_d(N_1+1:N))]/tau_d;
    
    % Dynamics
    xdd = 1/params.I_OS * ((-params.mass_OS*params.g*params.lc_OS*cos(x))+ FT_ext.*ma_ext + FT_flex.*ma_flex + act - 0.042*xd); %  + Tdamp + FT*ma);
    
    % backward euler
    % opti.subject_to(xd(1:N-1)*dt +x(1:N-1) == x(2:N));
    opti.subject_to(xd(1:N_1-1)*dt1 +x(1:N_1-1) == x(2:N_1));
    opti.subject_to(xd(N_1:N-1)*dt2 +x(N_1:N-1) == x(N_1+1:N));
    % opti.subject_to(xdd(1:N-1)*dt +xd(1:N-1) == xd(2:N));
    opti.subject_to(xdd(1:N_1-1)*dt1 +xd(1:N_1-1) == xd(2:N_1));
    opti.subject_to(xdd(N_1:N-1)*dt2 +xd(N_1:N-1) == xd(N_1+1:N));
    % opti.subject_to(dlMdt_ext(1:N-1)*dt + lMtilda_ext(1:N-1) == lMtilda_ext(2:N)); % VmTilde met factor 10 (MRS)
    opti.subject_to(dlMdt_ext(1:N_1-1)*dt1 + lMtilda_ext(1:N_1-1) == lMtilda_ext(2:N_1));
    opti.subject_to(dlMdt_ext(N_1:N-1)*dt2 + lMtilda_ext(N_1:N-1) == lMtilda_ext(N_1+1:N));
    % opti.subject_to(dlMdt_flex(1:N-1)*dt + lMtilda_flex(1:N-1) == lMtilda_flex(2:N)); % VmTilde met factor 10 (MRS)
    opti.subject_to(dlMdt_flex(1:N_1-1)*dt1 + lMtilda_flex(1:N_1-1) == lMtilda_flex(2:N_1));
    opti.subject_to(dlMdt_flex(N_1:N-1)*dt2 + lMtilda_flex(N_1:N-1) == lMtilda_flex(N_1+1:N));
    opti.subject_to(Fsrs_ddt(1:N_1-1)*dt1 + Fsrs_d(1:N_1-1) == Fsrs_d(2:N_1));
    opti.subject_to(Fsrs_ddt(N_1:N-1)*dt2 + Fsrs_d(N_1:N-1) == Fsrs_d(N_1+1:N));
    opti.subject_to(lM_ext.^2 - w_ext.^2 == lM_projected_ext.^2);
    opti.subject_to(lM_flex.^2 - w_flex.^2 == lM_projected_flex.^2);
    opti.subject_to(err_ext == 0);
    opti.subject_to(err_flex == 0);
    
    % Objective function
    tvect_spline = MX(1,N);
    for i = 1:N_1
        t = tvect(1)+dt1*i;
        tvect_spline(i) = t;
    end
    
    for i = N_1+1:N
        t = tvect(1)+dt1*N_1+dt2*i;
        tvect_spline(i) = t;
    end
    
    error        = x - q_exp;
    error_dot    = xd - qdot_exp;
    
    J            = sumsqr(error)  + sumsqr(error_dot) + 100*sumsqr(act); % + ...
    % 0.001 * (sumsqr(vMtilda_ext) + sumsqr(vMtilda_ext)); %+ 100*sumsqr(error_ra);
    opti.minimize(J);
    
    % options for IPOPT
    options.ipopt.tol = 1*10^(-6);          % moet normaal 10^-6 zijn
    options.ipopt.linear_solver = 'mumps';
    %     options.ipopt.nlp_scaling_method = 'none';
    %     options.ipopt.linear_solver = 'ma57';
    %options.ipopt.hessian_approximation = 'limited-memory'; % enkel bij moeilijke problemen
    
    % Solve the OCP
    opti.solver('ipopt',options);
    sol = opti.solve();
    
    sol_x = sol.value(x);
        sol_a_ext = sol.value(a_ext);               sol_a_flex = sol.value(a_flex);
    sol_lMtilda_ext = sol.value(lMtilda_ext);   sol_lMtilda_flex = sol.value(lMtilda_flex);
    sol_aext0 = sol.value(a_ext_0);
    sol_act = sol.value(act);
    sol_FT_ext  = sol.value(FT_ext);            sol_FT_flex  = sol.value(FT_flex);
    sol_ma_ext = sol.value(ma_ext);             sol_ma_flex = sol.value(ma_flex);
    sol_Fpe_ext = sol.value(Fpe_ext);           sol_Fpe_flex = sol.value(Fpe_flex);
    %sol_kFpe    = sol.value(kFpe);
    sol_kFpe_ext = sol.value(kFpe_ext);
    sol_kFpe_flex = sol.value(kFpe_flex); 
    %sol_kFpe_ext = sol.value(kFpe_ext);         sol_kFpe_flex = sol.value(kFpe_flex);
    sol_J   = sol.value(J);
    sol_Fsrs1 = sol.value(Fsrs1);
    sol_Fsrs2 = sol.value(Fsrs2);
    sol_Fsrs = [sol_Fsrs1 sol_Fsrs2];
    sol_Fsrs_d = sol.value(Fsrs_d);
    sol_dt1  = sol.value(dt1);
    sol_Rk   = sol.value(Rk);
    sol_dlmdt_ext = sol.value(dlMdt_ext);       sol_dlmdt_flex = sol.value(dlMdt_flex);
    sol_lM_ext = sol.value(lM_ext);             sol_lM_flex    = sol.value(lM_flex);
    sol_FMltilda_ext = sol.value(FMltilda_ext); sol_FMltilda_flex = sol.value(FMltilda_flex);
    Total = sol.value(FT_ext).*sol.value(ma_ext) + sol.value(FT_flex).*sol.value(ma_flex);
    save(['C:\Users\u0125183\Box\PhD 1\Simulations Pendulum Test\Results/Result_',char(s.nu),'_T',num2str(s.tr(j)),char(opt),'.mat'],'sol_aext0','sol_Rk','sol_dt1','sol_act','sol_x', 'sol_a_flex','sol_a_ext', 'sol_lMtilda_ext','sol_lMtilda_flex', 'sol_FT_ext', 'sol_FT_flex', 'sol_ma_ext', 'sol_ma_flex', 'sol_Fpe_ext', 'sol_Fpe_flex', 'sol_kFpe_ext','sol_kFpe_flex', 'sol_J', 'sol_Fsrs','q_exp','sol_dlmdt_ext','sol_dlmdt_flex','sol_lM_ext','sol_lM_flex','sol_FMltilda_ext','sol_FMltilda_flex','tvect')
    
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
end
