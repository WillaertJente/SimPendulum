%% Calculate ID based on results 
% To compare with experimental data
% Jente Willaert - 9 sept 2022

% Change with trial
load('CP11_T2_Opt14_Activation_TanH_tresh_FS_IG8.mat')
id_exp  = importdata(['C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Implicit\Muscle\Experimental data\New\CP11\ID_T2/inverse_dynamics.sto']);
emg_exp = xlsread(['C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\EMG/CP11/T2_Sim.xlsx']);

% input
m       = R.OS.inert.mass_OS;
g       = 9.81;
l       = R.OS.inert.lc_OS;
theta   = R.x;
FT_ext  = R.C.FT_ext;
FT_flex = R.C.FT_flex;
ma_ext  = R.C.MA_ext;
ma_flex = R.C.MA_flex;
B       = R.B;
theta_d = R.xd;

% ID
ID = FT_ext.*ma_ext + FT_flex.*ma_flex - B*theta_d; 
%-m*g*l*cos(theta*180/pi) + FT_ext.*ma_ext + FT_flex.*ma_flex - B*theta_d; 

% Spline experimental ID
tvect_spline = R.exp.t(1):R.dt2:R.exp.t(end);     % Time vector for simulation
N_spline     = length(tvect_spline) ;                      % Number of steps in simulation

idspline              = spline(id_exp.data(:,1),id_exp.data(:,R.subject.leg));            % Spline fit
[id_spline,qd_spline] = SplineEval_ppuval(idspline,tvect_spline,1); % Get angles and velocities

% Spline experimental emg
tvect_spline = R.exp.t(1):R.dt2:R.exp.t(end);     % Time vector for simulation
N_spline     = length(tvect_spline) ;                      % Number of steps in simulation

emgspline              = spline(emg_exp(:,1),emg_exp(:,2));            % Spline fit
[emg_spline,qd_spline] = SplineEval_ppuval(emgspline,tvect_spline,1); % Get angles and velocities

% Simulated emg
emg_sim = R.C.a;

% Plot
figure()
subplot(421)
plot(R.exp.qspline*180/pi,'k','LineWidth',1.5); hold on
plot(R.x*180/pi,'r','LineWidth',1.5); hold on
box off; title('CP11 T2'); 

subplot(4,2,3)
plot(id_spline,'k','LineWidth',1.5); hold on
plot(ID,'r','LineWidth',1.5); hold on
box off; title('ID')

subplot(4,2,4)
plot(FT_ext.*ma_ext,'LineWidth',1.5,'Color',[0 180 216]./255); hold on
plot(FT_flex.*ma_flex,'LineWidth',1.5,'Color',[179 146 172]./255); hold on
box off; title('Muscle moments'); legend({'ext','flex'})

subplot(4,2,5,'replace')
plot(emg_spline,'k','LineWidth',1.5); hold on
box off; title('experimental EMG')
subplot(4,2,7)
plot(emg_sim,'r','LineWidth',1.5); hold on
box off; title('Simulated EMG')

subplot(3,2,6,'replace')
bar(1,R.a_ext,'FaceColor',[0 180 216]./255,'EdgeColor',[0 180 216]./255); hold on
bar(2,R.a_flex,'FaceColor',[179 146 172]./255,'EdgeColor',[179 146 172]./255); hold on
bar(3,R.kR,'FaceColor',[0.7 0.7 0.7],'EdgeColor',[0.7 0.7 0.7]); hold on
box off; title('Neural params'); xticks([1 2 3]); xticklabels({'a Ext','aFlex','R'})

