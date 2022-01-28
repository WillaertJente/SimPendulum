%% Plot Different results 

map = 'C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/';
color = [0 0 0; 0 95 115; 10 147 150; 148 210 189; 238 155 0; 202 103 2; 187 62 3; 155 34 38]./255
figure()
for i = 2:8
    Name = [map,'TD5_T2_Opt7_Reflexes_Tau0.0',num2str(i),'.mat']; 
    load(Name)
    if i == 2
        plot(R.exp.qspline,'k','LineWidth',1.5); hold on
        plot(R.x,'Color',color(i,:)); hold on
    else
        plot(R.x,'Color',color(i,:)); hold on
    end
    title('TD5 - T1'); box off; 
    
end