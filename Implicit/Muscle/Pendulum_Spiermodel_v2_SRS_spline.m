%% Simulations of pendulum test using a muscle model
%  Jente Willaert - 21/10/2020
% clear all; close all; clc;

%% Input
s.nu = 'TD5';                                                                               % subject number/ name
s.tr = [2];                                                                                 % subject trials (number of trials)
pathmain = pwd;
[pathTemp,~,~] = fileparts(pathmain);
[pathRepo,~,~] = fileparts(pathTemp);
path = [pathRepo '\Implicit\Muscle\Experimental data\' s.nu '\'];                       % Path to opensim model (scaled)
ScaleFactor =1.697; % TD5 = 1.697 CP 4 = 1.5757 CP 8 = 1.7036
opt  = 'Spline2';   % Option used as name to save results

params = ImportParameters(s.nu);    % Input parameters (mtot,lc, l, age, m, RG, SE, Nmr, z)

for j = 1:length(s.tr)
    
    %% Prepare experimental data for simulation
    % Load experimental data + SRS on or off
    % Q_exp = BK data - 90°
    [q_exp_r, t_exp, t_span, on_srs] = LoadExpData(s.nu,s.tr,j,0, params, path);      % 1 if you want to plot experimental data
    
    % Discretised time: interpolate experimental data at discr. time using spline
    [q_exp, qdot_exp,N, tvect, dt] = ExpDatAtDiscrTime(t_span,t_exp,q_exp_r);       % dt = 0.005
    
%     q_exp = q_exp(1:1600);
%     qdot_exp = qdot_exp(1:1600);
%     N = 1600; 
%     tvect = tvect(1:1600);
%     tspan  = [tvect(1) tvect(end)]
    %% Formulate OCP
    import casadi.*;        % Import casadi libraries
    opti = casadi.Opti();   % Initialise opti structure
    xgrid = tvect;
    V     = q_exp;
    W     = qdot_exp;
    q_exp_spline    = casadi.interpolant('LUT','bspline',{xgrid},V);
    q_dot_spline    = casadi.interpolant('LUT','bspline',{xgrid},W);
    
    % Define phases of pendulum (initial state, end of first swing)
    [x0, N_1] = PendulumPhases(q_exp, qdot_exp, N, 1);                          % 1 if you want to plot the phases
    
    %% Add OpenSim Model
    addpath('MuscleModel');
    
    % Opensim model
    import org.opensim.modeling.*
    %model_path = 'C:\Opensim 3.3\Models\Gait2392_Simbody\gait2392_simbody.osim';
    model_path = [path,'/TD5_ScaledModel_ScaledForces.osim'];                                 % if cp = CPModel_Scaled.osim
    osimModel  = Model(model_path);
    
    % Inertial parameters (tibia)
    bodies         = osimModel.getBodySet();
    if params.z == 18
        tibia      = bodies.get('tibia_l');
    else
        tibia      = bodies.get('tibia_r');
    end
    params.mass_OS = tibia.getMass();
    
    com            = ArrayDouble.createVec3(0);
    tibia.getMassCenter(com);
    params.lc_OS   = abs(com.get(1));
    
    inertia        = Mat33(0);
    tibia.getInertia(inertia);
    params.I_OS    = inertia.get(0,0) + params.mass_OS * params.lc_OS ^2; %(apply steiner)
    
    % Muscle tendon properties
    [params_Muscle,lOpt,L_TendonSlack,Fiso,PennationAngle]=ReadMuscleParameters(model_path,{'rect_fem_r'});
    params_Muscle(2,1) = params_Muscle(2,1);
    params_Muscle(3,1) = params_Muscle(3,1);
    params.MTparams    = params_Muscle; % 5 x 2 matrix: 1 - Fmo, 2 - Lmo, 3 - Lts, 4 - alphao, 5 - vmmax
    
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
    
    % lM projected guess
    lMo      = params.MTparams(2,:)';
    alphao   = params.MTparams(4,:)';
    lMGuess  = lMtildeGuess.*lMo;
    w        = lMo.*sin(alphao);
    lM_projectedGuess = sqrt((lMGuess.^2 - w.^2));
    
    %dLdMT guess
    dlMdtGuess = vMGuess*params.MTparams(5)/lMo;
    
    %% Calculate LMT en Ma
    map_MA         = [path 'MA_FakeMot_T',num2str(s.tr(j))];
    [coeff_LMT_ma] = DefineLMTCoefficients(map_MA, s.nu);
    
    %% Calculate offset (difference between IK en BK)
    % BK + offset = IK
    IKtrial  = [path,'IK_Trial0',num2str(s.tr),'.mot'];
    [offset] = CalculateOffset(IKtrial);
    
    %% Define states, controls, bounds and initial guess
    
    % States
    x       = opti.variable(1,N);
    xd      = opti.variable(1,N);
    lMtilda = opti.variable(1,N);
    Fsrs    = opti.variable(1,N);
    
    % Controls
    lM_projected = opti.variable(1,N);
    act          = opti.variable(1,N);
    
    % Slack controls
    vMtilda = opti.variable(1,N);
    dt1          = opti.variable(1);
%     dt2          = opti.variable(1);
    
    % Parameters
    a            = opti.variable(1);
    kFpe         = opti.variable(1);
    
    % Bounds
    opti.subject_to(-4*pi < x   < 4*pi);
    opti.subject_to(-300  < xd  < 300);
    opti.subject_to(0.001 < a   < 0.5);
    opti.subject_to(1e-4  < lM_projected);       % Only positive lM's
    opti.subject_to(-10   < vMtilda < 10);
    opti.subject_to(0.2   < lMtilda < 1.8);
    opti.subject_to(0     < Fsrs    < 2);
    opti.subject_to(0.05  < kFpe    < 0.15);
    opti.subject_to(0.001 < dt1     < 0.01);    % 0.05
  %  opti.subject_to(0.001 < dt2     < 0.01);
    % opti.subject_to(act <= 0);
    
    % Bounds on initial states
    opti.subject_to(x(1)     == x0(1));
    opti.subject_to(xd(1)    == 0);
    %opti.subject_to(xd(1)    == x0(2));
    
    % Initial guess
    opti.set_initial(x, q_exp);
    opti.set_initial(xd,qdot_exp);
    opti.set_initial(a, 0.01);
    opti.set_initial(kFpe,0.1);
    opti.set_initial(lM_projected, lM_projectedGuess);
    opti.set_initial(lMtilda, lMtildeGuess);         % LmTildeGuess
    opti.set_initial(vMtilda, vMGuess);
    opti.set_initial(act,0);
    opti.set_initial(dt1,0.005);
    %.set_initial(dt2,0.005);
    
    % Defining problem (muscle model)
    % Calculate shift
    kT = 35;
    [shift]    = getshift(kT);
    
    % Calculate FT en ma
    [FT, ma, dlMdt, err, lM, lT,Fce, Fpe, FM, Fsrs, Fsrs_dot, FMltilda] = CalculateTendonForceAndMomentArm_v2_SRS(x, params, lMtilda, a, s.nu,shift, vMtilda, lM_projected,coeff_LMT_ma, offset, kFpe, N_1, Fsrs, N);
    
    % Variable time of phases
    tF1 = (N_1)*dt1;        %N_1-1
    dt2 = 0.005;
    
    % Constraints on phases
    opti.subject_to(tF1 > 0.2 );
    opti.subject_to(1e-4 < xd(N_1+1) < 0.1);
            
    % Constraints srs
    for k = 1:N_1
        opti.subject_to(xd(k) < 0);
    end
    for k = N_1:N-1
        opti.subject_to(Fsrs(k+1) == Fsrs(k) + Fsrs_dot(k) * dt2)
    end
    
    % Dynamics
    xdd = 1/params.I_OS * ((-params.mass_OS*params.g*params.lc_OS*cos(x))+ FT.*ma + act ); %  + Tdamp + FT*ma);
      
    % backward euler
    %     opti.subject_to(xd(1:N-1)*dt +x(1:N-1) == x(2:N));
    opti.subject_to(xd(1:N_1-1)*dt1 +x(1:N_1-1) == x(2:N_1));
    opti.subject_to(xd(N_1:N-1)*dt2 +x(N_1:N-1) == x(N_1+1:N));
    %    opti.subject_to(xdd(1:N-1)*dt +xd(1:N-1) == xd(2:N));
    opti.subject_to(xdd(1:N_1-1)*dt1 +xd(1:N_1-1) == xd(2:N_1));
    opti.subject_to(xdd(N_1:N-1)*dt2 +xd(N_1:N-1) == xd(N_1+1:N));
%    opti.subject_to(dlMdt(1:N-1)*dt + lMtilda(1:N-1) == lMtilda(2:N)); % VmTilde met factor 10 (MRS)
    opti.subject_to(dlMdt(1:N_1-1)*dt1 + lMtilda(1:N_1-1) == lMtilda(2:N_1));
    opti.subject_to(dlMdt(N_1:N-1)*dt2 + lMtilda(N_1:N-1) == lMtilda(N_1+1:N));
    opti.subject_to(lM.^2 - w.^2 == lM_projected.^2)
    opti.subject_to(err == 0);
    
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
       
     error        = x - q_exp_spline(tvect_spline);
     error_dot    = xd - q_dot_spline(tvect_spline);
    %  error = x-q_exp;
    %  error_dot = xd-qdot_exp;
    error_fs     = x(N_1)-q_exp(N_1);
    error_fs_dot = xd(N_1)-qdot_exp(N_1);
%    error_ra     = x(N-500:end)-q_exp(N-500:end);
    
    %error_fs = x(1:N_1+10) - q_exp(1:N_1+10);
    J        = sumsqr(error)  + sumsqr(error_dot) + sumsqr(act);
    opti.minimize(J);
    
    % options for IPOPT
    options.ipopt.tol = 1*10^(-6);          % moet normaal 10^-6 zijn
%         options.ipopt.linear_solver = 'mumps';
    options.ipopt.linear_solver = 'ma57';
    %options.ipopt.hessian_approximation = 'limited-memory'; % enkel bij moeilijke problemen
    
    % Solve the OCP
    opti.solver('ipopt',options);
    sol = opti.solve();
    
    sol_x = sol.value(x);
    sol_a = sol.value(a);
    sol_lMtilda = sol.value(lMtilda);
    sol_act = sol.value(act);
    sol_FT  = sol.value(FT);
    sol_ma = sol.value(ma);
    sol_Fpe = sol.value(Fpe);
    sol_kFpe = sol.value(kFpe);
    sol_J   = sol.value(J);
    sol_Fsrs = sol.value(Fsrs);
    save(['C:\Users\u0125183\Box\PhD 1\Simulations Pendulum Test\Results/Result_',char(s.nu),'_T',num2str(s.tr(j)),'_',char(opt),'.mat'],'sol_x', 'sol_a', 'sol_lMtilda', 'sol_act', 'sol_FT', 'sol_ma', 'sol_Fpe', 'sol_kFpe', 'sol_J', 'sol_Fsrs')
    
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
    plot(tvect,sol.value(FT),'LineWidth',1.5);
    hold on; box off; ylabel('FT');
    subplot(613)
    plot(tvect,sol.value(ma),'LineWidth',1.5);
    hold on; box off; ylabel('ma');
    subplot(614)
    plot(tvect,sol.value(Fpe),'LineWidth',1.5);
    hold on; box off; ylabel('Fpe');
    subplot(615)
    plot(tvect,sol.value(Fsrs),'LineWidth',1.5);
    hold on; box off; ylabel('Fsrs');
    subplot(616)
    plot(tvect,sol.value(act),'LineWidth',1.5);
    hold on; box off; ylabel('act');    
    
end
