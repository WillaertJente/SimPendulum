%% Plot differences between B values
%% CP 2
figure()
for i = 4:7
    name = ['T',num2str(i)];
    Name1 = ['C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results\CP2_T', num2str(i),'_Opt7_ingekort_Bopt.mat']
    Name2 = ['C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results\CP2_T', num2str(i),'_Opt7_Tau0.04.mat']
    CP2.Bh.(char(name)) = load(Name1);
    CP2.Bl.(char(name)) = load(Name2);
    
    subplot(2,2,i-3)
    plot(CP2.Bh.(char(name)).R.exp.qspline,'k','LineWidth',1.5); hold on
    plot(CP2.Bh.(char(name)).R.x,'r','LineWidth',1.5); hold on
    plot(CP2.Bl.(char(name)).R.x,'b','LineWidth',1.5); hold on
    legend({'Exp','Tau 0.08','Tau 0.04'}); box off; title(name)
end
suptitle('CP2')
saveas(gcf,'C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/Figures/CP2_InfluenceTau.fig')
saveas(gcf,'C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/Figures/CP2_InfluenceTau.png')

%% TD5
figure()
for i = 1:2
    name = ['T',num2str(i)];
    Name1 = ['C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results\TD5_T', num2str(i),'_Opt7_Bopt.mat']
    Name2 = ['C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results\TD5_T', num2str(i),'_Opt7_Bscaled.mat']
    TD5.Bh.(char(name)) = load(Name1);
    TD5.Bl.(char(name)) = load(Name2);
    
    subplot(2,2,i)
    plot(TD5.Bh.(char(name)).R.exp.qspline,'k','LineWidth',1.5); hold on
    plot(TD5.Bh.(char(name)).R.x,'r','LineWidth',1.5); hold on
    plot(TD5.Bl.(char(name)).R.x,'b','LineWidth',1.5); hold on
    legend({'Exp','B opt','B scaled'}); box off; title(name)
end
suptitle('TD5')
saveas(gcf,'C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/Figures/TD5_InfluenceBscaled.fig')
saveas(gcf,'C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/Figures/TD5_InfluenceBscaled.png')

%% CP8
figure()
for i = [2 3 4 5]
    name = ['T',num2str(i)];
    Name1 = ['C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results\CP8_T', num2str(i),'_Opt7_ingekort_Bopt.mat']
    Name2 = ['C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results\CP8_T', num2str(i),'_Opt7_Tau0.04.mat']
    CP8.Bh.(char(name)) = load(Name1);
    CP8.Bl.(char(name)) = load(Name2);
    
    subplot(2,2,i-1)
    plot(CP8.Bh.(char(name)).R.exp.qspline,'k','LineWidth',1.5); hold on
    plot(CP8.Bh.(char(name)).R.x,'r','LineWidth',1.5); hold on
    plot(CP8.Bl.(char(name)).R.x,'b','LineWidth',1.5); hold on
    legend({'Exp','Tau 0.08','Tau 0.04'}); box off; title(name)
end
suptitle('CP8')
saveas(gcf,'C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/Figures/CP8_InfluenceTau.fig')
saveas(gcf,'C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/Figures/CP8_InfluenceTau.png')

%% CP4
figure()
for i = 3:4
    name = ['T',num2str(i)];
    Name1 = ['C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results\CP4_T', num2str(i),'_Opt7_Bopt.mat']
    Name2 = ['C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results\CP4_T', num2str(i),'_Opt7_Bscaled.mat']
    CP4.Bh.(char(name)) = load(Name1);
    CP4.Bl.(char(name)) = load(Name2);
    
    subplot(2,2,i-2)
    plot(CP4.Bh.(char(name)).R.exp.qspline,'k','LineWidth',1.5); hold on
    plot(CP4.Bh.(char(name)).R.x,'r','LineWidth',1.5); hold on
    plot(CP4.Bl.(char(name)).R.x,'b','LineWidth',1.5); hold on
    legend({'Exp','B opt','B scaled'}); box off; title(name)
end
suptitle('CP4')
saveas(gcf,'C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/Figures/CP4_InfluenceBscaled.fig')
saveas(gcf,'C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/Figures/CP4_InfluenceBscaled.png')

%% CP9
figure()
for i = 3:5
    name = ['T',num2str(i)];
    Name1 = ['C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results\CP9_T', num2str(i),'_Opt7_Bopt.mat']
    Name2 = ['C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results\CP9_T', num2str(i),'_Opt7_Bscaled.mat']
    CP9.Bh.(char(name)) = load(Name1);
    CP9.Bl.(char(name)) = load(Name2);
    
    subplot(2,3,i-1)
    plot(CP9.Bh.(char(name)).R.exp.qspline,'k','LineWidth',1.5); hold on
    plot(CP9.Bh.(char(name)).R.x,'r','LineWidth',1.5); hold on
    plot(CP9.Bl.(char(name)).R.x,'b','LineWidth',1.5); hold on
    legend({'Exp','B opt','B scaled'}); box off; title(name)
end
suptitle('CP9')
saveas(gcf,'C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/Figures/CP9_InfluenceBscaled.fig')
saveas(gcf,'C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/Figures/CP9_InfluenceBscaled.png')

%% CP10
figure()
for i = 13:14
    name = ['T',num2str(i)];
    Name1 = ['C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results\CP10_T', num2str(i),'_Opt7_ingekort_Bopt.mat']
    Name2 = ['C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results\CP10_T', num2str(i),'_Opt7_Tau0.04.mat']
    CP10.Bh.(char(name)) = load(Name1);
    CP10.Bl.(char(name)) = load(Name2);
    
    subplot(2,2,i-12)
    plot(CP10.Bh.(char(name)).R.exp.qspline,'k','LineWidth',1.5); hold on
    plot(CP10.Bh.(char(name)).R.x,'r','LineWidth',1.5); hold on
    plot(CP10.Bl.(char(name)).R.x,'b','LineWidth',1.5); hold on
    legend({'Exp','B opt','B scaled'}); box off; title(name)
end
suptitle('CP10')
saveas(gcf,'C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/Figures/CP10_InfluenceTau.fig')
saveas(gcf,'C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/Figures/CP10_InfluenceTau.png')

%% CP11
figure()
for i = 2:4
    name = ['T',num2str(i)];
    Name1 = ['C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results\CP11_T', num2str(i),'_Opt7_ingekort_Bopt.mat']
    Name2 = ['C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results\CP11_T', num2str(i),'_Opt7_Tau0.04.mat']
    CP11.Bh.(char(name)) = load(Name1);
    CP11.Bl.(char(name)) = load(Name2);
    
    subplot(2,2,i-1)
    plot(CP11.Bh.(char(name)).R.exp.qspline,'k','LineWidth',1.5); hold on
    plot(CP11.Bh.(char(name)).R.x,'r','LineWidth',1.5); hold on
    plot(CP11.Bl.(char(name)).R.x,'b','LineWidth',1.5); hold on
    legend({'Exp','B opt','B scaled'}); box off; title(name)
end
suptitle('CP11')
saveas(gcf,'C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/Figures/CP11_InfluenceTau.fig')
saveas(gcf,'C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/Figures/CP11_InfluenceTau.png')