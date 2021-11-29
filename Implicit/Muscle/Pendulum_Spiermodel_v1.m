%% Simulations of pendulum test using a muscle model 
%  Jente Willaert - 21/10/2020
clear all; close all; clc; 

%% Input 
s.nu = 'TD6';                                                                               % subject number/ name 
s.tr = [1];                                                                                 % subject trials (number of trials)
ScaleFactor = 1.697 ; % TD5 = 1.697 CP 4 = 1.5757 CP 8 = 1.7036
opt  = '';                                                                                  % Option used as name to save results                                                                                % Weigths in cost function

pathmain = pwd;
[pathTemp,~,~] = fileparts(pathmain);
[pathRepo,~,~] = fileparts(pathTemp);
path = [pathRepo '\Implicit\Muscle\Experimental data\' s.nu '\'];                        % Path to opensim model (scaled)

params = ImportParameters(s.nu);    % Input parameters (mtot,lc, l, age, m, RG, SE, Nmr, z)

for j = 1:length(s.tr)
    
    %% Prepare experimental data for simulation
    % Load experimental data + SRS on or off
    % Q_exp = BK data - 90°
    [q_exp_r, t_exp, t_span, on_srs] = LoadExpData(s.nu,s.tr,j,0, params, path);      % 1 if you want to plot experimental data
    
    % Discretised time: interpolate experimental data at discr. time using spline 
    [q_exp, qdot_exp,N, tvect, dt] = ExpDatAtDiscrTime(t_span,t_exp,q_exp_r);       % dt = 0.005

    % Define phases of pendulum (initial state, end of first swing) 
    [x0, N_1] = PendulumPhases(q_exp, qdot_exp, N, 1);                          % 1 if you want to plot the phases

    %% Formulate OCP 
    import casadi.*;        % Import casadi libraries
    opti = casadi.Opti();   % Initialise opti structure
    
    %% Add OpenSim Model
    addpath('MuscleModel');
    
    % Opensim model
    import org.opensim.modeling.*           
    %model_path = 'C:\Opensim 3.3\Models\Gait2392_Simbody\gait2392_simbody.osim'; 
    model_path = [path,'TD5_ScaledModel_ScaledForces.osim'];                                 % if cp = CPModel_Scaled.osim
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
    xdd     = opti.variable(1,N);
    lMtilda = opti.variable(1,N); 
    dlMdt   = opti.variable(1,N);
    vMtilda = opti.variable(1,N);
        
    % Controls
    a            = opti.variable(1);
    lM_projected = opti.variable(1,N);
    act          = opti.variable(1,N);
    kFpe         = opti.variable(1); 
%     klim         = opti.variable(1); 
%     Kr1_OS       = opti.variable(1); 
%     Kr2_OS       = opti.variable(1); 

%     % Bounds
    opti.subject_to(-4*pi < x   < 4*pi);  
    opti.subject_to(-300  < xd  < 300);
    opti.subject_to(0.001 < a   < 1);
    opti.subject_to(1e-4  < lM_projected);       % Only positive lM's
    opti.subject_to(-10   < vMtilda < 10);
    opti.subject_to(0.2   < lMtilda < 1.5);
    opti.subject_to(-0.2  < kFpe    < 0.2); 
%     opti.subject_to(-160  < Kr1_OS < -120);
%     opti.subject_to(-30   < Kr2_OS < 10);
%     opti.subject_to(0     < klim    < 6);
    opti.subject_to(act <= 0);

%     
    % Bounds on initial states
    opti.subject_to(x(1)     == x0(1));
    %opti.subject_to(xd(1)    == x0(2));
    opti.subject_to(xd(1)    == 0);
    
    % Initial guess
    opti.set_initial(x, q_exp);
    opti.set_initial(xd,qdot_exp);
    opti.set_initial(a, 0.01);
    opti.set_initial(kFpe,0.1);
    opti.set_initial(lM_projected, lM_projectedGuess);
    opti.set_initial(lMtilda, lMtildeGuess);         % LmTildeGuess
    opti.set_initial(vMtilda, vMGuess); 
    opti.set_initial(dlMdt,dlMdtGuess);
    opti.set_initial(act,0); 
%     opti.set_initial(Kr1_OS,-140); 
%     opti.set_initial(Kr2_OS,-10); 
%     opti.set_initial(klim,10);
%        
    % Defining problem (muscle model)
    % Calculate shift
    kT = 35; 
    [shift]    = getshift(kT);                   
%     
    % Calculate Tlim
    klim   = 5;                                     % 3 in script torques, 10 in script Friedl 
    Kr1_OS = -140;
    Kr2_OS = -10;
%     TLim = 0;
   %[TLim] = CalculateTLim(x,Kr1_OS, Kr2_OS, klim);
   [TLim] = CalculateTLim_KA(x,Kr1_OS, Kr2_OS, klim, offset);
   TLim = TLim/ScaleFactor;
    %[Tlim] = CalculateTLim(x, Kr1, Kr2, klim);
    
    % Calculate FT en ma 
    [FT, ma, dlMdt, err, lM, lT,Fce, Fpe, FM] = CalculateTendonForceAndMomentArm(x, params, lMtilda, a, s.nu,shift, vMtilda, lM_projected,coeff_LMT_ma, offset, kFpe);
    
    
    % Dynamics
    xdd = 1/params.I_OS * ((-params.mass_OS*params.g*params.lc_OS*cos(x))+ FT.*ma  + TLim + act); %  + Tdamp + FT*ma);
    
    % backward euler
    opti.subject_to(xd(1:N-1)*dt +x(1:N-1) == x(2:N));
    opti.subject_to(xdd(1:N-1)*dt +xd(1:N-1) == xd(2:N));
    opti.subject_to(dlMdt(1:N-1)*dt + lMtilda(1:N-1) == lMtilda(2:N)); % VmTilde met factor 10 (MRS)
    opti.subject_to(lM.^2 - w.^2 == lM_projected.^2)
    opti.subject_to(err == 0); 
%     
    % Objective function
    error = x - q_exp; 
    J     = sumsqr(error) + 2*sumsqr(act);
    opti.minimize(J); 
    
    
    % options for IPOPT
    options.ipopt.tol = 1*10^(-10);          % moet normaal 10^-6 zijn
    options.ipopt.linear_solver = 'mumps';
%     options.ipopt.linear_solver = 'ma57';
    %options.ipopt.hessian_approximation = 'limited-memory'; % enkel bij moeilijke problemen
          
    % Solve the OCP
    opti.solver('ipopt',options);
    sol = opti.solve();  
    
    sol_x = sol.value(x);
    
    figure(j*10)
    plot(tvect,q_exp,'k','LineWidth',1.5)
    hold on
    plot(tvect,sol_x,'r','LineWidth',1.5)
    hold on
    ylabel('Knee angle (rad)'); 
    hold on
    title('Check tracking');
    xlabel('Time (s)');
    legend('Experimental data','Simulated data')
    
    figure(j*100)
    subplot(711)
    plot(tvect,q_exp,'k','LineWidth',1.5)
    hold on
    plot(tvect,sol_x,'LineWidth',1.5)
    hold on
    ylabel('Knee angle (rad)'); box off; legend('Exp','Sim');
    subplot(712)
    plot(tvect,sol.value(TLim),'LineWidth',1.5);
    hold on; box off; ylabel('TLim (N)');
    subplot(713)
    plot(tvect,sol.value(FT),'LineWidth',1.5);
    hold on; box off; ylabel('FT (N)');         % Tendon force
    subplot(714)
    plot(tvect,sol.value(ma),'LineWidth',1.5);  % Moment arm
    hold on; box off; ylabel('ma (m)');
    subplot(715)
    plot(tvect,sol.value(Fpe),'LineWidth',1.5); % Passive force
    hold on; box off; ylabel('Fpe (N)');
    subplot(716)
    scatter(tvect(50),sol.value(kFpe),'LineWidth',1.5);
    hold on; box off; ylabel('kFpe');           % Optimized passive force constant
    subplot(717)
    plot(tvect,sol.value(act),'LineWidth',1.5);
    hold on; box off; ylabel('actuator (Nm)'); xlabel('Time (s)')

end
