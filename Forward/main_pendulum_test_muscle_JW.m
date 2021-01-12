% Simulate a pendulum test.
%
% Friedl De Groote
% August 24, 2017

% Jente Willaert UPDATA
% August 2020
% Code works until the part without SRS and with feedback. 
% Feedback is work in progress 

% clear all
% close all

addpath('MuscleModel');

% Colors
color_palette= [255 0 0; 255 150 150; 0 114 54; 153 255 102; 0 0 255; 44 176 207; 0 0 0; 155 155 155; 255 153 0; 255 185 0; 153 51 102; 171 218 77]*(1/255);
color=color_palette([6 1 2 10 11 4 3 5 12 7 8 9 1 10 4 6 11 3 2 5 12],:); %S0052 functional vars (EMGGRFandCOM file)

%% Get mass and inertia parameters from OpenSim model
% Import the OpenSim modeling classes
import org.opensim.modeling.*

% Let's use an OpenSim model
% For now, a simple generic model
model_path = 'C:\OpenSim 3.3\Models\Gait2392_Simbody/gait2392_simbody.osim';
%'C:\OpenSim 3.3\Models\Gait10dof18musc/gait10dof18musc.osim';


% Read in the osim model
osimModel = Model(model_path);

% Get the inertial parameters
% For now, let's assume that the tibia is most important
% methodsview(osimModel);
bodies = osimModel.getBodySet();
tibia = bodies.get('tibia_r');

% Get the mass
inputdata.mass = tibia.getMass();

% Get center of mass location
com = ArrayDouble.createVec3(0);
tibia.getMassCenter(com);
inputdata.lc = abs(com.get(1));

% Get the inertia
inertia = Mat33(0);
tibia.getInertia(inertia);
inputdata.I = inertia.get(0,0) + inputdata.mass * inputdata.lc ^ 2; % Apply Steiner

% Read in muscle-tendon properties
% [params,lOpt,L_TendonSlack,Fiso,PennationAngle]=ReadMuscleParameters(model_path,{'rect_fem_r', 'vasti_r'});
[params,lOpt,L_TendonSlack,Fiso,PennationAngle]=ReadMuscleParameters(model_path,{'rect_fem_r'});
params(2,1) = params(2,1);
params(3,1) = params(3,1);
inputdata.MTparams = params; % 5 x 2 matrix: 1 - Fmo, 2 - Lmo, 3 - Lts, 4 - alphao, 5 - vmmax

% Load muscle-tendon lengths
LMTdata = importdata('LMT_knee_ext_92.mot');
inputdata.LMT_qknee = [-160; -140; LMTdata.data(:,2)] * pi/180;
% inputdata.LMT = [LMTdata.data(1,3:4); LMTdata.data(1,3:4); LMTdata.data(:,3:4)]; % RF - VAS
inputdata.LMT = [LMTdata.data(1,3); LMTdata.data(1,3); LMTdata.data(:,3)]; % RF

% Load moment arms
madata = importdata('ma_knee_ext_92.mot');
inputdata.ma_qknee = [-160; -140; madata.data(:,2)] * pi/180;
% inputdata.ma = [madata.data(1,3:4); madata.data(1,3:4); madata.data(:,3:4)]; % RF - VAS
inputdata.ma = [madata.data(1,3); madata.data(1,3); madata.data(:,3)]; % RF - VAS

% Parameters of active muscle force-velocity characteristic
load Fvparam.mat
Fvparam(1) = 1.475*Fvparam(1); Fvparam(2) = 0.25*Fvparam(2);
Fvparam(3) = Fvparam(3) + 0.75; Fvparam(4) = Fvparam(4) - 0.027;
inputdata.Fvparam = Fvparam;

% Parameters of active muscle force-length characteristic
load Faparam.mat                            
inputdata.Faparam = Faparam;

% Parameters of passive muscle force-length characteristic
e0 = 0.6; kpe = 4; t50 = exp(kpe * (0.2 - 0.10e1) / e0);
pp1 = (t50 - 0.10e1); t7 = exp(kpe); pp2 = (t7 - 0.10e1);
inputdata.Fpparam = [pp1;pp2];

%% Simulate pendulum movement
global hist
% parameters
inputdata.d = 0.175; %0.350  % 0.175
inputdata.knee_r_range = [-120 0]*pi/180;  
inputdata.klim = 10; %10
inputdata.kM = 3;   % Maxwell 
inputdata.kF = 10;
inputdata.kdF = 0.02;
inputdata.kl = 0.02; % 
inputdata.kv = 0; % 0.01; 
% inputdata.Tb = 1;

% SRS
% inputdata.delta_theta_crit = 5/180 * pi;
% inputdata.kSRS = 1/inputdata.delta_theta_crit;
inputdata.ksrs = 1000;
inputdata.Tb = 2;
inputdata.tauSRS = 0.05;

tspan = [0 12]';
kin0 = [0/180*pi 0]';
x0srs = [0/180*pi 0 0 0]';
options = [];
lags = 0.05;
inputdata.lags = lags;


%% Simulate pendulum movement (torque driven)
% Maxwell model
inputdata.kM = 3;   % Maxwell
inputdata.d  = 10;

inputdata.Tb = 0.3;
[tM,xM] = ode15s(@pendulumStateDerivativeMKV, tspan, [kin0; 0], options, inputdata);

figure()
plot(tM,xM(:,1)*180/pi, 'Color', [0 0 0], 'LineWidth',2); hold on;
xlabel('time [s]')
ylabel('angle [^o]')

legend_names = {'passive joint torque'};
i = 2;

%% Simulate pendulum movement WITH srs (torque driven) (Jente)

% Inputdata damping
inputdata.kM = 3;   % Maxwell
inputdata.d  = 10;

% Inputdata Tb
inputdata.Tb = 0.3;

% inputdata voor Tsrs
inputdata.start = kin0(1);
inputdata.qc    = -1.5*pi/180;         % Qc = 1.5 graad = critical angle (Tsrs = ct*delta q if stretch is smaller than Qc / Tsrs = ct * Qc if stretch is larger than Qc)
inputdata.kSRS  = 0.67*180/pi; 
inputdata.tauSRS = 0.05;

hist = 1
[tM,xM] = ode15s(@pendulumStateDerivativeMKVsrs, tspan, [kin0; 0; 0], options, inputdata);

figure()
plot(tM,xM(:,1)*180/pi, 'Color', [0 0 0], 'LineWidth',2); hold on;
xlabel('time [s]')
ylabel('angle [^o]')

legend_names = {'passive joint torque'};
i = 2;

%% Simulate pendulum movement (muscle driven)
% inputdata.ab = 0.01*ones(1,2); % baseline muscle tone
inputdata.ab = 0.01; %0.01
inputdata.kM = 3;   % Maxwell 
inputdata.d  = 0.175;

% Generate initial guess of lMtilda
inputdata.q0 = kin0(1);
[tiso,xiso] = ode15s(@isometricDerivative, [0 0.5], [0.3], options, inputdata);
xiso(end)
x0 = [kin0; 0; xiso(end,:)'];

% Simulate pendulum movement
[tM,xM] = ode15s(@pendulumStateDerivativeMuscle, tspan, x0, options, inputdata);

plot(tM,xM(:,1)*180/pi, 'Color', color(i-1,:), 'LineWidth',2); hold on;
legend_names{i} = 'muscle - a_b = 0.01';
i = i+1;

% Increased baseline tone
% inputdata.ab = 0.03*ones(1,2); % baseline muscle tone
inputdata.ab = 0.02;
% Generate initial guess of lMtilda
inputdata.q0 = kin0(1);
[tiso,xiso] = ode15s(@isometricDerivative, [0 0.5], [0.3], options, inputdata);
xiso(end,:)
x0 = [kin0; 0; xiso(end,:)'];

% Simulate pendulum movement
[tM,xM] = ode15s(@pendulumStateDerivativeMuscle, tspan, x0, options, inputdata);

plot(tM,xM(:,1)*180/pi, 'Color', color(i-1,:), 'LineWidth',2); hold on;
legend_names{i} = 'muscle - a_b = 0.02';
i = i+1;

%% Simulate pendulum movement (muscle driven) with SRS (Jente aangepast)
% inputdata.ab = 0.01*ones(1,2); % baseline muscle tone
inputdata.ab   = 0.06;
inputdata.d    = 0.175;
inputdata.ksrs = 150;
inputdata.tauSRS = 0.05;
inputdata.kM   = 3; 

% Generate initial guess of lMtilda
inputdata.q0 = kin0(1);
[tiso,xiso] = ode15s(@isometricDerivative, [0 0.5],0.3, options, inputdata);  %0.3
inputdata.lMiso = xiso(end,:)'; 
x0 = [kin0; 0; xiso(end,:)';0];  % laatste 0 toegevoegd? 

% Simulate pendulum movement
hist = 1;
[tM,xM] = ode15s(@pendulumStateDerivativeMuscleSRS_JW, tspan, x0, options, inputdata);

plot(tM,xM(:,1)*180/pi, 'Color', color(i-1,:), 'LineWidth',2); hold on;
legend_names{i} = 'muscle SRS - a_b = 0.02';
i = i+1

%% Simulate pendulum movement (muscle driven) withOUT SRS and delayed feedback
% % inputdata.ab = 0.01*ones(1,2); % baseline muscle tone
% inputdata.ab = 0.1;
% inputdata.d = 0.175;
% inputdata.klim = 1000;
% lags = 0.1;
% inputdata.lags = lags;
% 
% % Generate initial guess of lMtilda
% inputdata.q0 = kin0(1);
% [tiso,xiso] = ode15s(@isometricDerivative, [0 0.5], [0.3], options, inputdata);
% 
% % Compute normalized isometric force
% LMT = interp1(inputdata.LMT_qknee, inputdata.LMT, inputdata.q0);
% [dx,FT] = FiberLengthOde(inputdata.ab,xiso(end,:)',LMT, inputdata.MTparams, inputdata.Fvparam, inputdata.Fpparam, inputdata.Faparam);
% inputdata.fse_iso = FT./inputdata.MTparams(1,:);
% inputdata.lMiso = xiso(end,:)'; 
% x0 = [kin0; 0; inputdata.ab; xiso(end,:)'];
% 
% % Simulate pendulum movement
% hist = 0;
% global params
% params = inputdata;
% solFsSRS = ddensd(@pendulumStateDerivativeMuscleSRS_forceFB, lags, lags, x0, tspan);
% 
% plot(solFsSRS.x,solFsSRS.y(1,:)*180/pi, 'Color',color(i-1,:),'LineWidth',2); hold on;
% iFsSRS = i;
% legend_names{i} = 'muscle - force FB - a_b = 0.02';
% i = i+1

%% Simulate pendulum movement (muscle driven) withOUT SRS and delayed feedback (JENTE)
% inputdata.ab = 0.01*ones(1,2); % baseline muscle tone
inputdata.ab = 0.06;
inputdata.d = 0.175;
inputdata.klim = 10;
lags = 0.1;
inputdata.lags = lags;

% Generate initial guess of lMtilda
inputdata.q0 = kin0(1);
[tiso,xiso] = ode15s(@isometricDerivative, [0 0.5], 0.3, options, inputdata);

% Compute normalized isometric force
LMT = interp1(inputdata.LMT_qknee, inputdata.LMT, inputdata.q0);
[dx,FT] = FiberLengthOde(inputdata.ab,xiso(end,:)',LMT, inputdata.MTparams, inputdata.Fvparam, inputdata.Fpparam, inputdata.Faparam);
inputdata.fse_iso = FT./inputdata.MTparams(1,:);
inputdata.lMiso = xiso(end,:)'; 
x0 = [kin0; 0; inputdata.ab; xiso(end,:)'];

% Simulate pendulum movement
hist = 0;
global params
params = inputdata;
solFsSRS = ddensd(@pendulumStateDerivativeMuscleSRS_forceFB, lags, lags, x0, tspan);

plot(solFsSRS.x,solFsSRS.y(1,:)*180/pi, 'Color',color(i-1,:),'LineWidth',2); hold on;
iFsSRS = i;
legend_names{i} = 'muscle - force FB - a_b = 0.02';
i = i+1

%% Simulate pendulum movement (muscle driven) with SRS and delayed feedback
% inputdata.ab = 0.01*ones(1,2); % baseline muscle tone
inputdata.ab = 0.01;
inputdata.d  = 0.175;
inputdata.klim = 10;
inputdata.ksrs = 1000;
inputdata.tauSRS = 0.05;
% Generate initial guess of lMtilda
inputdata.q0 = kin0(1);
[tiso,xiso] = ode15s(@isometricDerivative, [0 0.5], [0.3], options, inputdata);

% Compute normalized isometric force
LMT = interp1(inputdata.LMT_qknee, inputdata.LMT, inputdata.q0);
[dx,FT] = FiberLengthOde(inputdata.ab,xiso(end,:)',LMT, inputdata.MTparams, inputdata.Fvparam, inputdata.Fpparam, inputdata.Faparam);
inputdata.fse_iso = FT./inputdata.MTparams(1,:);
inputdata.lMiso = xiso(end,:)'; 
x0 = [kin0; 0; inputdata.ab; xiso(end,:)'];

% Simulate pendulum movement
hist = 1;
global params
params = inputdata;
solFsSRS = ddensd(@pendulumStateDerivativeMuscleSRS_forceFB, lags, lags, x0, tspan);

plot(solFsSRS.x,solFsSRS.y(1,:)*180/pi, 'Color',color(i-1,:),'LineWidth',2); hold on;
iFsSRS = i;
legend_names{i} = 'muscle SRS - force FB - a_b = 0.02';
i = i+1

%% Simulate pendulum movement (muscle driven) with SRS and delayed feedback
% inputdata.ab = 0.01*ones(1,2); % baseline muscle tone
inputdata.ab = 0.02;
% Generate initial guess of lMtilda
inputdata.q0 = kin0(1);
[tiso,xiso] = ode15s(@isometricDerivative, [0 0.5], [0.3], options, inputdata);
inputdata.lMiso = xiso(end,:)'; 
x0 = [kin0; 0; inputdata.ab; xiso(end,:)'];

% Simulate pendulum movement
hist = 0;
global params
params = inputdata;
sollsSRS = ddensd(@pendulumStateDerivativeMuscleSRS_lengthFB, lags, lags, x0, tspan);

plot(sollsSRS.x,sollsSRS.y(1,:)*180/pi, 'Color',color(i-1,:),'LineWidth',2); hold on;
ilsSRS = i;
legend_names{i} = 'muscle - length FB - a_b = 0.02';
i = i+1;

legend(legend_names)

%% Simulate pendulum movement (muscle driven) with SRS and delayed feedback
% inputdata.ab = 0.01*ones(1,2); % baseline muscle tone
inputdata.ab = 0.02;
% Generate initial guess of lMtilda
inputdata.q0 = kin0(1);
[tiso,xiso] = ode15s(@isometricDerivative, [0 0.5], [0.3], options, inputdata);
inputdata.lMiso = xiso(end,:)'; 
x0 = [kin0; 0; inputdata.ab; xiso(end,:)'];

% Simulate pendulum movement
hist = 1;
global params
params = inputdata;
sollsSRS = ddensd(@pendulumStateDerivativeMuscleSRS_lengthFB, lags, lags, x0, tspan);

plot(sollsSRS.x,sollsSRS.y(1,:)*180/pi, 'Color',color(i-1,:),'LineWidth',2); hold on;
ilsSRS = i;
legend_names{i} = 'muscle SRS - length FB - a_b = 0.02';
i = i+1;

legend(legend_names)

% plot muscle activation
figure()
plot(solFsSRS.x,solFsSRS.y(4,:), 'Color',color(iFsSRS-1,:),'LineWidth',2); hold on;
plot(sollsSRS.x,sollsSRS.y(4,:), 'Color',color(ilsSRS-1,:),'LineWidth',2); hold on;
