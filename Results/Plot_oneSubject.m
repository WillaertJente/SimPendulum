%% CP 2
figure()
for i = 4:7
    name = ['T',num2str(i)];
    Name = ['CP2_T', num2str(i),'_Opt7_BInitValue.mat']
    CP2.(char(name)) = load(Name);
    
    subplot(2,5,i-3)
    plot(CP2.(char(name)).R.exp.qspline,'k','LineWidth',1.5); hold on
    plot(CP2.(char(name)).R.x,'r','LineWidth',1.5); hold on
    title(name); box off;
    ylim([-2 0]); xlim([0 1000])
    
    subplot(2,5,6)
    bar(i-3,CP2.(char(name)).R.a_ext); hold on
    title('a_ext'); box off; xticks([1 2 3 4]); xticklabels({'T4','T5','T6','T7'})
    ylim([0 0.02]);
    
    subplot(2,5,7)
    bar(i-3,CP2.(char(name)).R.a_flex); hold on
    title('a_flex'); box off; xticks([1 2 3 4]); xticklabels({'T4','T5','T6','T7'})
    ylim([0 0.1]);
    
    subplot(2,5,8)
    bar(i-3,CP2.(char(name)).R.kFpe); hold on
    title('kFpe'); box off; xticks([1 2 3 4]); xticklabels({'T4','T5','T6','T7'})
    ylim([0 0.2]); line([0 5],[0.2 0.2],'Color',[0.7 0.7 0.7],'LineWidth',1.5)
    
    subplot(2,5,9)
    bar(i-3,CP2.(char(name)).R.kR); hold on
    title('kR'); box off; xticks([1 2 3 4]); xticklabels({'T4','T5','T6','T7'})
    ylim([0 10]); line([0 5],[10 10],'Color',[0.7 0.7 0.7],'LineWidth',1.5)
    
    subplot(2,5,10)
    bar(i-3,CP2.(char(name)).R.B); hold on
    title('B'); box off; xticks([1 2 3 4]); xticklabels({'T4','T5','T6','T7'})
    ylim([0 0.12]); line([0 5],[10 10],'Color',[0.7 0.7 0.7],'LineWidth',1.5)
end
suptitle('CP2')
saveas(gcf,'C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/Figures/CP2_BInitValue.fig')
saveas(gcf,'C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/Figures/CP2_BInitValue.png')

%% CP 4
figure()
for i = 3:4
    name = ['T',num2str(i)];
    Name = ['CP4_T', num2str(i),'_Opt7_BInitValue.mat']
    CP4.(char(name)) = load(Name);
    
    subplot(2,5,i-2)
    plot(CP4.(char(name)).R.exp.qspline,'k','LineWidth',1.5); hold on
    plot(CP4.(char(name)).R.x,'r','LineWidth',1.5); hold on
    title(name); box off;
    ylim([-2.5 0]); xlim([0 2000])
    
    subplot(2,5,6)
    bar(i-2,CP4.(char(name)).R.a_ext); hold on
    title('a_ext'); box off; xticks([1 2]); xticklabels({'T3','T4'})
    ylim([0 0.02]);
    
    subplot(2,5,7)
    bar(i-2,CP4.(char(name)).R.a_flex); hold on
    title('a_flex'); box off; xticks([1 2]); xticklabels({'T3','T4'})
    ylim([0 0.1])
    
    subplot(2,5,8)
    bar(i-2,CP4.(char(name)).R.kFpe); hold on
    title('kFpe'); box off; xticks([1 2]); xticklabels({'T3','T4'})
    ylim([0 0.2]); line([0 3],[0.2 0.2],'Color',[0.7 0.7 0.7],'LineWidth',1.5)
    
    subplot(2,5,9)
    bar(i-2,CP4.(char(name)).R.kR); hold on
    title('kR'); box off; xticks([1 2]); xticklabels({'T3','T4'})
    ylim([0 10]); line([0 3],[10 10],'Color',[0.7 0.7 0.7],'LineWidth',1.5)
    
    subplot(2,5,10)
    bar(i-2,CP4.(char(name)).R.B); hold on
    title('B'); box off; xticks([1 2]); xticklabels({'T3','T4'})
    ylim([0 0.12]); line([0 3],[10 10],'Color',[0.7 0.7 0.7],'LineWidth',1.5)
end
suptitle('CP4')
saveas(gcf,'C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/Figures/CP4_BInitValue.mat.fig')
saveas(gcf,'C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/Figures/CP4_BInitValue.mat.png')

%% CP8
figure()
for i = [2 3 4 5]
    name = ['T',num2str(i)];
    Name = ['CP8_T', num2str(i),'_Opt7_BInitValue.mat']
    CP8.(char(name)) = load(Name);
    
    subplot(2,5,i-1)
    plot(CP8.(char(name)).R.exp.qspline,'k','LineWidth',1.5); hold on
    plot(CP8.(char(name)).R.x,'r','LineWidth',1.5); hold on
    title(name); box off;
    ylim([-2 0]); xlim([0 1500])
    
    subplot(2,5,6)
    bar(i-1,CP8.(char(name)).R.a_ext); hold on
    title('a_ext'); box off; xticks([1 2 3 4]); xticklabels({'T2','T3','T4','T5'})
    ylim([0 0.02]);
    
    subplot(2,5,7)
    bar(i-1,CP8.(char(name)).R.a_flex); hold on
    title('a_flex'); box off; xticks([1 2 3 4]); xticklabels({'T2','T3','T4','T5'})
    ylim([0 0.1])
    
    subplot(2,5,8)
    bar(i-1,CP8.(char(name)).R.kFpe); hold on
    title('kFpe'); box off; xticks([1 2 3 4]); xticklabels({'T2','T3','T4','T5'})
    ylim([0 0.2]); line([0 5],[0.2 0.2],'Color',[0.7 0.7 0.7],'LineWidth',1.5)
    
    subplot(2,5,9)
    bar(i-1,CP8.(char(name)).R.kR); hold on
    title('kR'); box off; xticks([1 2 3 4]); xticklabels({'T2','T3','T4','T5'})
    ylim([0 10]); line([0 5],[10 10],'Color',[0.7 0.7 0.7],'LineWidth',1.5)
    
    subplot(2,5,10)
    bar(i-1,CP8.(char(name)).R.B); hold on
    title('B'); box off; xticks([1 2 3 4]); xticklabels({'T2','T3','T4','T5'})
    ylim([0 0.12]); line([0 5],[10 10],'Color',[0.7 0.7 0.7],'LineWidth',1.5)
end
suptitle('CP8')
saveas(gcf,'C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/Figures/CP8_BInitValue.fig')
saveas(gcf,'C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/Figures/CP8_BInitValue.png')

%% CP9
figure()
for i = 2:6
    name = ['T',num2str(i)];
    Name = ['CP9_T', num2str(i),'_Opt7_BInitValue.mat']
    CP9.(char(name)) = load(Name);
    
    subplot(2,5,i-1)
    plot(CP9.(char(name)).R.exp.qspline,'k','LineWidth',1.5); hold on
    plot(CP9.(char(name)).R.x,'r','LineWidth',1.5); hold on
    title(name); box off;
    ylim([-2.5 0]); xlim([0 2000])
    
    subplot(2,5,6)
    bar(i-1,CP9.(char(name)).R.a_ext); hold on
    title('a_ext'); box off; xticks([1 2 3 4 5]); xticklabels({'T3','T4','T5','T6'})
    ylim([0 0.02]);
    
    subplot(2,5,7)
    bar(i-1,CP9.(char(name)).R.a_flex); hold on
    title('a_flex'); box off; xticks([1 2 3 4 5]); xticklabels({'T3','T4','T5','T6'})
    ylim([0 0.1])
    
    subplot(2,5,8)
    bar(i-1,CP9.(char(name)).R.kFpe); hold on
    title('kFpe'); box off; xticks([1 2 3 4 5]); xticklabels({'T3','T4','T5','T6'})
    ylim([0 0.2]); line([0 5],[0.2 0.2],'Color',[0.7 0.7 0.7],'LineWidth',1.5)
    
    subplot(2,5,9)
    bar(i-1,CP9.(char(name)).R.kR); hold on
    title('kR'); box off; xticks([1 2 3 4 5]); xticklabels({'T3','T4','T5','T6'})
    ylim([0 10]); line([0 5],[10 10],'Color',[0.7 0.7 0.7],'LineWidth',1.5)
    
    subplot(2,5,10)
    bar(i-1,CP9.(char(name)).R.B); hold on
    title('B'); box off; xticks([1 2 3 4 5]); xticklabels({'T3','T4','T5','T6'})
    ylim([0 0.12]); line([0 5],[10 10],'Color',[0.7 0.7 0.7],'LineWidth',1.5)
end
suptitle('CP9')
saveas(gcf,'C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/Figures/CP9_BInitValue.fig')
saveas(gcf,'C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/Figures/CP9_BInitValue.png')

%% CP 10
figure()
for i = 12:14
    name = ['T',num2str(i)];
    Name = ['CP10_T', num2str(i),'_Opt7_BInitValue.mat']
    CP10.(char(name)) = load(Name);
    
    subplot(2,5,i-11)
    plot(CP10.(char(name)).R.exp.qspline,'k','LineWidth',1.5); hold on
    plot(CP10.(char(name)).R.x,'r','LineWidth',1.5); hold on
    title(name); box off;
    ylim([-2 0]); xlim([0 1200])
    
    subplot(2,5,6)
    bar(i-11,CP10.(char(name)).R.a_ext); hold on
    title('a_ext'); box off; xticks([1 2 3]); xticklabels({'T12','T13','T14'})
    ylim([0 0.05]);
    
    subplot(2,5,7)
    bar(i-11,CP10.(char(name)).R.a_flex); hold on
    title('a_flex'); box off; xticks([1 2 3]); xticklabels({'T12','T13','T14'})
    ylim([0 0.1])
    
    subplot(2,5,8)
    bar(i-11,CP10.(char(name)).R.kFpe); hold on
    title('kFpe'); box off; xticks([1 2 3]); xticklabels({'T12','T13','T14'})
    ylim([0 0.2]); line([0 4],[0.2 0.2],'Color',[0.7 0.7 0.7],'LineWidth',1.5)
    
    subplot(2,5,9)
    bar(i-11,CP10.(char(name)).R.kR); hold on
    title('kR'); box off; xticks([1 2 3]); xticklabels({'T12','T13','T14'})
    ylim([0 10]); line([0 4],[10 10],'Color',[0.7 0.7 0.7],'LineWidth',1.5)
    
    subplot(2,5,10)
    bar(i-11,CP10.(char(name)).R.B); hold on
    title('B'); box off; xticks([1 2 3]); xticklabels({'T12','T13','T14'})
    ylim([0 0.12]); line([0 4],[10 10],'Color',[0.7 0.7 0.7],'LineWidth',1.5)
end
suptitle('CP10')
saveas(gcf,'C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/Figures/CP10_BInitValue.fig')
saveas(gcf,'C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/Figures/CP10_BInitValue.png')

%% CP 11
figure()
for i = 2:4
    name = ['T',num2str(i)];
    Name = ['CP11_T', num2str(i),'_Opt7_BInitValue.mat']
    CP11.(char(name)) = load(Name);
    
    subplot(2,5,i-1)
    plot(CP11.(char(name)).R.exp.qspline,'k','LineWidth',1.5); hold on
    plot(CP11.(char(name)).R.x,'r','LineWidth',1.5); hold on
    title(name); box off;
    ylim([-2.5 0]); xlim([0 1500])
    
    subplot(2,5,6)
    bar(i-1,CP11.(char(name)).R.a_ext); hold on
    title('a_ext'); box off; xticks([1 2]); xticklabels({'T3','T4'})
    ylim([0 0.02]);
    
    subplot(2,5,7)
    bar(i-1,CP11.(char(name)).R.a_flex); hold on
    title('a_flex'); box off; xticks([1 2]); xticklabels({'T3','T4'})
    ylim([0 0.1])
    
    subplot(2,5,8)
    bar(i-1,CP11.(char(name)).R.kFpe); hold on
    title('kFpe'); box off; xticks([1 2]); xticklabels({'T3','T4'})
    ylim([0 0.2]); line([0 3],[0.2 0.2],'Color',[0.7 0.7 0.7],'LineWidth',1.5)
    
    subplot(2,5,9)
    bar(i-1,CP11.(char(name)).R.kR); hold on
    title('kR'); box off; xticks([1 2]); xticklabels({'T3','T4'})
    ylim([0 10]); line([0 3],[10 10],'Color',[0.7 0.7 0.7],'LineWidth',1.5)
    
    subplot(2,5,10)
    bar(i-1,CP11.(char(name)).R.B); hold on
    title('B'); box off; xticks([1 2]); xticklabels({'T3','T4'})
    ylim([0 0.12]); line([0 3],[10 10],'Color',[0.7 0.7 0.7],'LineWidth',1.5)
end
suptitle('CP11')
saveas(gcf,'C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/Figures/CP11_BInitValue.fig')
saveas(gcf,'C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/Figures/CP11_BInitValue.png')

%% TD5
figure()
for i = 1:2
    name = ['T',num2str(i)];
    Name = ['TD5_T', num2str(i),'_Opt7_BInitValue.mat']
    TD5.(char(name)) = load(Name);
    
    subplot(2,5,i)
    plot(TD5.(char(name)).R.exp.qspline,'k','LineWidth',1.5); hold on
    plot(TD5.(char(name)).R.x,'r','LineWidth',1.5); hold on
    title(name); box off;
    ylim([-2.5 0]); xlim([0 2500])
    
    subplot(2,5,6)
    bar(i,TD5.(char(name)).R.a_ext); hold on
    title('a_ext'); box off; xticks([1 2]); xticklabels({'T1','T2'})
    ylim([0 0.02]);
    
    subplot(2,5,7)
    bar(i,TD5.(char(name)).R.a_flex); hold on
    title('a_flex'); box off; xticks([1 2]); xticklabels({'T1','T2'})
    ylim([0 0.1])
    
    subplot(2,5,8)
    bar(i,TD5.(char(name)).R.kFpe); hold on
    title('kFpe'); box off; xticks([1 2]); xticklabels({'T1','T2'})
    ylim([0 0.2]); line([0 3],[0.2 0.2],'Color',[0.7 0.7 0.7],'LineWidth',1.5)
    
    subplot(2,5,9)
    bar(i,TD5.(char(name)).R.kR); hold on
    title('kR'); box off; xticks([1 2]); xticklabels({'T1','T2'})
    ylim([0 0.1]); line([0 3],[10 10],'Color',[0.7 0.7 0.7],'LineWidth',1.5)
    
    subplot(2,5,10)
    bar(i,TD5.(char(name)).R.B); hold on
    title('B'); box off; xticks([1 2]); xticklabels({'T1','T2'})
    ylim([0 0.12]); line([0 3],[10 10],'Color',[0.7 0.7 0.7],'LineWidth',1.5)
end
suptitle('TD5')
saveas(gcf,'C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/Figures/TD5_BInitValue.mat.fig')
saveas(gcf,'C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/Figures/TD5_BInitValue.mat.png')
