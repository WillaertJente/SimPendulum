function [h] = PlotResults(R, q_forward, info)
% Plot results of forward and tracking simulations

% Info on figure
lw = 1.5; % linewidth
h  = figure('Name','Simulation results');

plot(R.exp.qspline*180/pi,'k','LineWidth',lw); hold on
plot(R.x*180/pi,'r','LineWidth',lw); hold on
plot(q_forward*180/pi,'Color',[35 87  137]./255, 'LineWidth', lw); hold on

legend({'Experimental','Tracking','Forward'},'Location','Best');
ylabel('Knee angle ({\circ})')
xlabel('Frames (200Hz)')
title([info.subj ' Trial ', num2str(info.trial)])
box off; 
end