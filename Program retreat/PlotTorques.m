function p = PlotTorques(R)
% Plot torques

% Info on figure
lw = 1.5; % linewidth
p  = figure('Name','Simulated torques');

subplot(411)
plot(R.x*180/pi,'k','LineWidth',lw); hold on
plot(R.xd*180/pi,'Color',[0.6 0.6 0.6],'LineWidth',lw); hold on
plot(R.xdd*180/pi,'Color',[0.8 0.8 0.8],'LineWidth',lw); hold on
legend({'x','xd','xdd'}); title('Simulated trajectory'); ylabel('[{\circ}]');box off;

subplot(412)
plot(R.C.T_Fz,'k','LineWidth',1.5); hold on
title('Torque Fz'); ylabel('[Nm]'); box off;

subplot(413)
plot(R.C.T_inert,'k','LineWidth',1.5); hold on
title('Torque Inertia'); ylabel('[Nm]'); box off; 

subplot(414)
plot(R.C.T_muscle,'k','LineWidth',1.5); hold on
title('Torque muscle'); ylabel('[Nm]'); box off 
end

