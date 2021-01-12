%% Plot knee angle, velocity and acceleration

% Importdata
map = 'C:\Users\u0125183\Documents\PhD 1\Accelerometer\Vicon/';

for i = [10]
% Plot knee angle 

    if i < 10
        d = importdata([map,'Trial0',num2str(i),'.mot']); 
    else
        d = importdata([map,'Trial',num2str(i),'.mot']); 
    end
    
    figure(i)
    subplot(211)
    plot(d.data(:,11)*pi/180,'LineWidth',1.5)
    hold on
    ylabel('Knee angle (rad)')
    title(['Trial',num2str(i)])
    
% Define drop based on velocity 
    v = diff(d.data(:,18)*pi/180)./(d.data(2,1)-d.data(1,1));
    subplot(212)
    plot(v,'LineWidth',1.5)
    hold on
    
    dropf = find(v <= -0.1);
    drop  = dropf(17);
    line([drop,drop],[-0.5 0.5],'Color',[0 0 0],'LineWidth',1.5)
    hold on
    ylabel('velocity (rad/s)')
    
    subplot(211)
    line([drop,drop],[-3 0],'Color',[0 0 0],'LineWidth',1.5)
    
end
