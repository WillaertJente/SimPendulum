function [h] = PlotResults(R, info)
% Plot results of forward and tracking simulations

% Info on figure
lw = 1.5; % linewidth
h  = figure('Name','Simulation results');

plot(R.exp.qspline*180/pi,'k','LineWidth',lw); hold on
plot(R.x*180/pi,'r','LineWidth',lw); hold on

legend({'Experimental','Tracking'},'Location','Best');
ylabel('Knee angle ({\circ})')
xlabel('Frames (200Hz)')
title([info.subj ' Trial ', num2str(info.trial)])
box off; 
end