%% Plot EMG with results optimization
color = [0 48 73; 214 40 40; 247 127 0; 252 191 73; 234 226 183]./255; 
k = 1
% Results = Opt7_Ingekort_BOpt
subj = 'CP9';
figure()
for t = 2:6
    map_res = ['C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/',char(subj),'_T',num2str(t),'_Opt7_ingekort_Bopt.mat'];
    map_emg = ['C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\EMG/',char(subj),'/T',num2str(t),'_Sim.xlsx'];
    
    res = load(map_res);
    emg = xlsread(map_emg,'Sheet1');
    
    subplot(3,5,k)
    plot(res.R.exp.tspline, res.R.exp.qspline,'k','LineWidth',1.5); hold on
    plot(res.R.exp.tspline,res.R.x,'r','LineWidth',1.5); hold on
    box off; title(['T',num2str(t)]); xlim([emg(1,1) emg(end,1)+1]); 
    
    subplot(3,5,5+k)
    plot(emg(:,1),emg(:,2),'k','LineWidth',1.5); hold on
    box off; xlim([emg(1,1) emg(end,1)+1]); ylim([0 0.00001])
    
    subplot(3,5,11)
    bar(k,res.R.a_ext,'FaceColor',color(k,:),'EdgeColor',color(k,:)); hold on
    box off; title('A Ext'); xticks(1:5); ylim([0 0.02])
    
    subplot(3,5,12)
    bar(k,res.R.a_flex,'FaceColor',color(k,:),'EdgeColor',color(k,:)); hold on
    box off; title('A Flex'); xticks(1:5); ylim([0 0.02])
    
    subplot(3,5,13)
    bar(k,res.R.kR,'FaceColor',color(k,:),'EdgeColor',color(k,:)); hold on
    box off; title('kR'); xticks(1:5); ylim([0 1])
    
    subplot(3,5,14)
    bar(k,res.R.kFpe,'FaceColor',color(k,:),'EdgeColor',color(k,:)); hold on
    box off; title('kFpe'); xticks(1:5); ylim([0 0.2])
    
    subplot(3,5,15)
    bar(k,res.R.B,'FaceColor',color(k,:),'EdgeColor',color(k,:)); hold on
    box off; title('B'); xticks(1:5); ylim([0 0.15])
    
    k = k+1
end

saveas(gcf,['C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results\Figures/',char(subj),'_EMG.png'])
saveas(gcf,['C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results\Figures/',char(subj),'_EMG.fig'])