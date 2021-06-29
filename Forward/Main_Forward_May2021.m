%% Forward implementation of pendulum test
% Adaptation of Friedl De Groote - August 24, 2017
% Jente Willaert may 27 2021
clear all; % close all;

vOpenSim = 4;

pathmain = pwd;
[pathTemp,~,~] = fileparts(pathmain);
[pathRepo,~,~] = fileparts(pathTemp);

% Path and subject information
subj     = 'TD5';
s.tr     = [1];
paramss  = ImportParameters_forward(subj);
pathmain = pwd;
addpath('MuscleModel');
% load('C:\Users\u0125183\Box\PhD 1\Simulations Pendulum Test\Results/Result_TD5_T1_Fmv.mat')
load('Result_TD5_T1_Fmv.mat')

% Colors
color_palette= [255 0 0; 255 150 150; 0 114 54; 153 255 102; 0 0 255; 44 176 207; 0 0 0; 155 155 155; 255 153 0; 255 185 0; 153 51 102; 171 218 77]*(1/255);
color=color_palette([6 1 2 10 11 4 3 5 12 7 8 9 1 10 4 6 11 3 2 5 12],:); %S0052 functional vars (EMGGRFandCOM file)

%% Load experimental trial
% Load experimental data
path = [pathRepo '\SimPendulum\Implicit\Muscle\Experimental data/',subj,'/'];
[q_exp_r, t_exp, t_span, on_srs] = LoadExpData_forward(subj,s.tr,0,paramss,path);      % 1 if you want to plot experimental data

% Discretised time: interpolate experimental data at discr. time using spline
[q_exp, qdot_exp,N, tvect, dt] = ExpDatAtDiscrTime_forward(t_span,t_exp,q_exp_r);

%% Get mass and inertia parameters from OpenSim model
% Import the OpenSim modeling classes
import org.opensim.modeling.*

% Subject-specific scaled model
model_path = [pathRepo '\SimPendulum\Implicit\Muscle\Experimental data/',subj,'/',subj,'_ScaledModel_ScaledForces.osim'];

% Read in the osim model
osimModel = Model(model_path);

if vOpenSim == 3
    % Get the inertial parameters (changed with before, now tibia + calcaneus)
    bodies         = osimModel.getBodySet();
    if paramss.z == 18
        tibia      = bodies.get('tibia_l');
        calcn      = bodies.get('calcn_l');
    else
        tibia      = bodies.get('tibia_r');
        calcn      = bodies.get('calcn_r');
    end
    
    % Get mass
    params.mass_OStibia = tibia.getMass();
    params.mass_OScalcn = calcn.getMass();
    params.mass_OS   = params.mass_OStibia + params.mass_OScalcn;
    
    % Get COM
    com_tibia         = ArrayDouble.createVec3(0);     tibia.getMassCenter(com_tibia);
    com_calcn         = ArrayDouble.createVec3(0);     calcn.getMassCenter(com_calcn);
    com_foot_in_tibia = abs(com_calcn.get(1)) + paramss.length_tibia;
    params.lc_OS   = (abs(com_tibia.get(1)) * tibia.getMass() + com_foot_in_tibia*calcn.getMass())/params.mass_OS;
    
    % Get inertia
    inertia_tibia  = Mat33(0);     tibia.getInertia(inertia_tibia);
    inertia_calcn  = Mat33(0);     calcn.getInertia(inertia_calcn);
    params.I_OS = inertia_tibia.get(0,0) + inertia_calcn.get(0,0) + tibia.getMass() * com_tibia.get(1) ^ 2 + calcn.getMass() * com_foot_in_tibia ^ 2;
    
    % Read in muscle-tendon properties
    [params_Muscle,lOpt,L_TendonSlack,Fiso,PennationAngle]=ReadMuscleParameters(model_path,{'rect_fem_r','bifemlh_r'});
    params_Muscle_ext(2,1) = params_Muscle(2,1);  params_Muscle_flex(2,1) = params_Muscle(2,2);
    params_Muscle_ext(3,1) = params_Muscle(3,1);  params_Muscle_flex(3,1) = params_Muscle(3,2);
    params.MTparams_ext    = [params_Muscle(:,1); sol_kFpe_ext];
    params.MTparams_flex   = [params_Muscle(:,2); sol_kFpe_flex];% 5 x 2 matrix: 1 - Fmo, 2 - Lmo, 3 - Lts, 4 - alphao, 5 - vmmax

elseif vOpenSim == 4
    % Get the inertial parameters (changed with before, now tibia + calcaneus)
    bodies         = osimModel.getBodySet();
    if paramss.z == 18
        tibia      = bodies.get('tibia_l');
        calcn      = bodies.get('calcn_l');
    else
        tibia      = bodies.get('tibia_r');
        calcn      = bodies.get('calcn_r');
    end
    
    % Get mass
    params.mass_OStibia = tibia.getMass();
    params.mass_OScalcn = calcn.getMass();
    params.mass_OS   = params.mass_OStibia + params.mass_OScalcn;
    
    % Get COM
    com_tibia         = tibia.getMassCenter();
    com_calcn         = calcn.getMassCenter();
    com_foot_in_tibia = abs(com_calcn.get(1)) + paramss.length_tibia;
    params.lc_OS   = (abs(com_tibia.get(1)) * tibia.getMass() + com_foot_in_tibia*calcn.getMass())/params.mass_OS;
    
    % Get inertia
    inertia_tibia  = tibia.getInertia();
    diagI_tibia = inertia_tibia.getMoments;
    offDiagI_tibia = inertia_tibia.getProducts;
    
    inertia_calcn  = calcn.getInertia();
    diagI_calcn = inertia_calcn.getMoments;
    offDiagI_calcn = inertia_calcn.getProducts;
    
    params.I_OS = diagI_tibia.get(0) + diagI_calcn.get(0) + tibia.getMass() * com_tibia.get(1) ^ 2 + calcn.getMass() * com_foot_in_tibia ^ 2;
    
    % Read in muscle-tendon properties
    [params_Muscle,lOpt,L_TendonSlack,Fiso,PennationAngle]=ReadMuscleParameters(model_path,{'rect_fem_r','bifemlh_r'});
    params_Muscle_ext(2,1) = params_Muscle(2,1);  params_Muscle_flex(2,1) = params_Muscle(2,2);
    params_Muscle_ext(3,1) = params_Muscle(3,1);  params_Muscle_flex(3,1) = params_Muscle(3,2);
    params.MTparams_ext    = [params_Muscle(:,1); sol_kFpe_ext];
    params.MTparams_flex   = [params_Muscle(:,2); sol_kFpe_flex];% 5 x 2 matrix: 1 - Fmo, 2 - Lmo, 3 - Lts, 4 - alphao, 5 - vmmax
    
end

% Load muscle-tendon lengths
map_MA         = [pathRepo '\SimPendulum\Implicit\Muscle\Experimental data/',subj,'/MA_FakeMot_T',num2str(s.tr)];
[coeff_LMT_ma_ext, coeff_LMT_ma_flex] = DefineLMTCoefficients_2muscles_forward(map_MA, subj,0); % Set last input to 1 if you want to visualize the fit.
params.coeff_LMT_ma_ext  = coeff_LMT_ma_ext;
params.coeff_LMT_ma_flex = coeff_LMT_ma_flex;

% Parameters of active muscle force-velocity characteristic
load Fvparam.mat
Fvparam(1) = 1.475*Fvparam(1); Fvparam(2) = 0.25*Fvparam(2);
Fvparam(3) = Fvparam(3) + 0.75; Fvparam(4) = Fvparam(4) - 0.027;
params.Fvparam = Fvparam;

% Parameters of active muscle force-length characteristic
load Faparam.mat
params.Faparam = Faparam;

% Parameters of passive muscle force-length characteristic
e0 = 0.6; kpe = 4; t50 = exp(kpe * (0.2 - 0.10e1) / e0);
pp1 = (t50 - 0.10e1); t7 = exp(kpe); pp2 = (t7 - 0.10e1);
params.Fpparam = [pp1;pp2];

%% Inputdata
% global hist
% parameters which we want to change (initial - we take parameters from
% implicit solution)
inputdata.B         = paramss.B;                       % Damping (personalized value)
inputdata.ab_ext    = 0.041; % sol_aext0;                       % Baseline muscle tone extensor
inputdata.ab_flex   = sol_a_flex;                      % Baseline muscle tone flexor
inputdata.Rk        = 0; %sol_Rk;                          % Reflex gain extensor
inputdata.kY        = 0;                             % reflex gain yank
inputdata.act       = sol_act;

% SRS
inputdata.ksrs   = 280;
inputdata.tauSRS = 0.05;
inputdata.tau_d  = 0.05;

% Initial state q en qdot
tspan = t_span';
kin0 = [q_exp(1) qdot_exp(1)]'; % moet snelheid nul zijn om te beginnen?

global hist

%% Simulate pendulum movement - muscle driven WITH srs
options = odeset('RelTol',1e-6,'AbsTol',1e-6);
% Initial guess of lMtilda voor flexor en extensor
inputdata.q0 = kin0(1);
% Friedl - there is a problem here. Fiber lengthts can become negative.
% Should look into this.
[~,xiso_ext]  = ode15s(@dldtIsometric, [0 1], 0.7, options,params, inputdata,'ext');
[~,xiso_flex] = ode15s(@dldtIsometric, [0 1], 0.7, options, params, inputdata, 'flex');
inputdata.lMiso_ext  = xiso_ext(end,:);
inputdata.lMiso_flex = xiso_flex(end,:);

% Initial state [q qdot lMtilde_ext Fsrs Fsrs_delayed lMtilde_flex Ysrs_delayed]
x0 = [kin0; inputdata.lMiso_ext';0; 0; inputdata.lMiso_flex';0];

% Simulate pendulum movement
hist = 1;
[tM,xM] = ode15s(@pendulumStateDerivativeMuscleSRS_Yfb, tspan, x0, options, inputdata, params);

% Plot results
figure()
subplot(311)
plot(tM,xM(:,1)*180/pi,'LineWidth',2); hold on; % q
plot(tvect,q_exp*180/pi,'k','LineWidth',1.5); hold on;
plot(tvect,sol_x*180/pi,'r','LineWidth',1.5); hold on;
title('Trajectory'); legend({'Forward','Experimental','Implicit'})
subplot(312)
plot(tM,xM(:,4),'LineWidth',1.5); hold on % Fsrs
plot(tM,xM(:,5),'LineWidth',1.5); hold on; % Fsrs_delayed
title('SRS'); legend({'Fsrs ext','Fsrs delayed'})

%% Calculate FT_ext and FT_flex
FT_ext = computeFTfromqandlM(params.MTparams_ext,params.coeff_LMT_ma_ext,xM(:,1),xM(:,3));
FT_flex = computeFTfromqandlM(params.MTparams_flex,params.coeff_LMT_ma_flex,xM(:,1),xM(:,6));

subplot(313)
plot(tM,FT_ext,'LineWidth',1.5); hold on % Tendon force extensor
plot(tM,FT_flex,'LineWidth',1.5); hold on % Tendon force flexor
title('Tendon force'); legend({'Ext','Flex'})

