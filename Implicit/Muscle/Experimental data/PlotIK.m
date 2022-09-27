%% Plot IK 

map = 'C:\Users\u0125183\OneDrive - KU Leuven\Pendulum test\Data\CP11\IK/'

%% ALL
for i = [2 3 4 5 6]
    if i < 10
        Name = [map,'Trial0',num2str(i),'.mot']; 
    else
        Name = [map,'Trial',num2str(i),'.mot']; 
    end
    if exist(Name)
        data = importdata(Name)
        
        figure()
        plot(data.data(:,11),'LineWidth',1.5)
        title(['Trial',num2str(i)])
    end
end


%% SIMULATED
for i = [2 3 4]
    if i < 10 
        Name = [map,'Trial0',num2str(i),'.mot']; 
    else
        Name = [map,'Trial',num2str(i),'.mot']; 
    end
    if exist(Name)
        data = importdata(Name)
        
        figure(100)
        subplot(121)
        plot(data.data(:,11),'LineWidth',1.5); hold on %11
        box off; legend({'T2','T3','T4'}); title('Knee angle')
        
        subplot(122)
        plot(data.data(:,8),'LineWidth',1.5); hold on %% 8
        box off; title('Hip flexion angle'); ylim([0 30])
    end
end