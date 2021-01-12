%% Controleren of BK hoek van tibia = IK hoek + offset 
% OFFset = 90 - HF - PT
map = 'C:\Users\u0125183\Box\PhD 1\Simulations Pendulum Test\TD5/';
%% Load BK
BK            = importdata([map, '/BK_Trial01__BodyKinematics_pos_global.sto']);
col_tibia     = find(strcmp(BK.colheaders,'tibia_r_Oz'));
theta_tib_hor = BK.data(:,col_tibia);

figure()
plot(BK.data(:,1),theta_tib_hor-90,'k','LineWidth',1.5)
hold on

%% Load IK 
IK           = importdata([map, 'IK_Trial05.mot']);
col_knee     = find(strcmp(IK.colheaders,'knee_angle_r'));
theta_knee   = IK.data(:,col_knee);

plot(IK.data(:,1),theta_knee,'Color',[0 0 0.7],'LineWidth',1.5)
hold on

%% Calculate Offset
col_pelT   = find(strcmp(IK.colheaders,'pelvis_tilt'));
col_hipF   = find(strcmp(IK.colheaders,'hip_flexion_r'));

start      = ones(length(IK.data),1)*90;
offset     = start - IK.data(:,col_pelT) - IK.data(:,col_hipF);

controle   = theta_knee - offset;

plot(IK.data(:,1), controle, 'r','LineWidth',1.5)
hold on;
plot(IK.data(:,1),offset,'Color',[0.7 0.7 0.7],'LineWidth',1.5);
hold on; box off; ylabel('{\circ}'); xlabel('Time(s)');
legend('theta-Hor','knee flex','kneeF + offset')
saveas(gcf, [map,'ControleBK_Trial05.png'])