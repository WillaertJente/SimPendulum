%% Plot IK 

map = 'C:\Users\u0125183\Box\PhD 1\Simulations Pendulum Test\Data\CP1/'

for i = 1:13
    if i < 10
        Name = [map,'IK_Trial0',num2str(i),'.mot']; 
    else
        Name = [map,'IK_Trial',num2str(i),'.mot']; 
    end
    if exist(Name)
        data = importdata(Name)
        
        figure()
        plot(data.data(:,11),'LineWidth',1.5)
    end
end
