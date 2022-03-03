%% CP 2
color4 = [0 48 73]./255;
color5 = [214 40 40]./255;
color6 = [247 127 0]./255;
color7 = [252 191 73]./255;

figure()
trials = [4 5; 4 6; 4 7;5 6; 5 7; 6 7];
for i = 1:length(trials)
    
    name  = ['T',num2str(trials(i,1)),'-T',num2str(trials(i,2))];
    Name  = ['CP2_', name,'_Opt8_MT.mat']
    fname = ['V',num2str(i)];
    CP2.(char(fname)) = load(Name);
    
    subplot(2,6,i)
    plot(CP2.(char(fname)).R.exp.qspline,'k','LineWidth',1.5); hold on
    plot(CP2.(char(fname)).R.x,'r','LineWidth',1.5); hold on
    title(name); box off;
    ylim([-2 0]); xlim([0 1500])
end

subplot(2,6,7)
bar(1,CP2.V1.R.a_ext(1),'FaceColor',color4,'EdgeColor',color4); hold on
bar(2,CP2.V1.R.a_ext(2),'FaceColor',color5,'EdgeColor',color5); hold on
bar(4,CP2.V2.R.a_ext(1),'FaceColor',color4,'EdgeColor',color4); hold on
bar(5,CP2.V2.R.a_ext(2),'FaceColor',color6,'EdgeColor',color6); hold on
bar(7,CP2.V3.R.a_ext(1),'FaceColor',color4,'EdgeColor',color4); hold on
bar(8,CP2.V3.R.a_ext(2),'FaceColor',color7,'EdgeColor',color7); hold on
bar(10,CP2.V4.R.a_ext(1),'FaceColor',color5,'EdgeColor',color5); hold on
bar(11,CP2.V4.R.a_ext(2),'FaceColor',color6,'EdgeColor',color6); hold on
bar(13,CP2.V5.R.a_ext(1),'FaceColor',color5,'EdgeColor',color5); hold on
bar(14,CP2.V5.R.a_ext(2),'FaceColor',color7,'EdgeColor',color7); hold on
bar(16,CP2.V6.R.a_ext(1),'FaceColor',color6,'EdgeColor',color6); hold on
bar(17,CP2.V6.R.a_ext(2),'FaceColor',color7,'EdgeColor',color7); hold on
xticks([1 2 4 5 7 8 10 11 13 14 16 17]); xticklabels({'4','5','4','6','4','7','5','6','5','7','6','7'});
box off; title('a ext'); ylim([0 0.006])

subplot(2,6,8)
bar(1,CP2.V1.R.a_flex(1),'FaceColor',color4,'EdgeColor',color4); hold on
bar(2,CP2.V1.R.a_flex(2),'FaceColor',color5,'EdgeColor',color5); hold on
bar(4,CP2.V2.R.a_flex(1),'FaceColor',color4,'EdgeColor',color4); hold on
bar(5,CP2.V2.R.a_flex(2),'FaceColor',color6,'EdgeColor',color6); hold on
bar(7,CP2.V3.R.a_flex(1),'FaceColor',color4,'EdgeColor',color4); hold on
bar(8,CP2.V3.R.a_flex(2),'FaceColor',color7,'EdgeColor',color7); hold on
bar(10,CP2.V4.R.a_flex(1),'FaceColor',color5,'EdgeColor',color5); hold on
bar(11,CP2.V4.R.a_flex(2),'FaceColor',color6,'EdgeColor',color6); hold on
bar(13,CP2.V5.R.a_flex(1),'FaceColor',color5,'EdgeColor',color5); hold on
bar(14,CP2.V5.R.a_flex(2),'FaceColor',color7,'EdgeColor',color7); hold on
bar(16,CP2.V6.R.a_flex(1),'FaceColor',color6,'EdgeColor',color6); hold on
bar(17,CP2.V6.R.a_flex(2),'FaceColor',color7,'EdgeColor',color7); hold on
xticks([1 2 4 5 7 8 10 11 13 14 16 17]); xticklabels({'4','5','4','6','4','7','5','6','5','7','6','7'});
box off; title('a flex'); ylim([0 0.006])

subplot(2,6,9)
bar(1,CP2.V1.R.kR(1),'FaceColor',color4,'EdgeColor',color4); hold on
bar(2,CP2.V1.R.kR(2),'FaceColor',color5,'EdgeColor',color5); hold on
bar(4,CP2.V2.R.kR(1),'FaceColor',color4,'EdgeColor',color4); hold on
bar(5,CP2.V2.R.kR(2),'FaceColor',color6,'EdgeColor',color6); hold on
bar(7,CP2.V3.R.kR(1),'FaceColor',color4,'EdgeColor',color4); hold on
bar(8,CP2.V3.R.kR(2),'FaceColor',color7,'EdgeColor',color7); hold on
bar(10,CP2.V4.R.kR(1),'FaceColor',color5,'EdgeColor',color5); hold on
bar(11,CP2.V4.R.kR(2),'FaceColor',color6,'EdgeColor',color6); hold on
bar(13,CP2.V5.R.kR(1),'FaceColor',color5,'EdgeColor',color5); hold on
bar(14,CP2.V5.R.kR(2),'FaceColor',color7,'EdgeColor',color7); hold on
bar(16,CP2.V6.R.kR(1),'FaceColor',color6,'EdgeColor',color6); hold on
bar(17,CP2.V6.R.kR(2),'FaceColor',color7,'EdgeColor',color7); hold on
xticks([1 2 4 5 7 8 10 11 13 14 16 17]); xticklabels({'4','5','4','6','4','7','5','6','5','7','6','7'});
box off; title('kR'); ylim([0 2])

subplot(2,6,10)
bar(1,CP2.V1.R.kFpe,'FaceColor',color4,'EdgeColor',color4); hold on
bar(2,CP2.V1.R.kFpe,'FaceColor',color5,'EdgeColor',color5); hold on
bar(4,CP2.V2.R.kFpe,'FaceColor',color4,'EdgeColor',color4); hold on
bar(5,CP2.V2.R.kFpe,'FaceColor',color6,'EdgeColor',color6); hold on
bar(7,CP2.V3.R.kFpe,'FaceColor',color4,'EdgeColor',color4); hold on
bar(8,CP2.V3.R.kFpe,'FaceColor',color7,'EdgeColor',color7); hold on
bar(10,CP2.V4.R.kFpe,'FaceColor',color5,'EdgeColor',color5); hold on
bar(11,CP2.V4.R.kFpe,'FaceColor',color6,'EdgeColor',color6); hold on
bar(13,CP2.V5.R.kFpe,'FaceColor',color5,'EdgeColor',color5); hold on
bar(14,CP2.V5.R.kFpe,'FaceColor',color7,'EdgeColor',color7); hold on
bar(16,CP2.V6.R.kFpe,'FaceColor',color6,'EdgeColor',color6); hold on
bar(17,CP2.V6.R.kFpe,'FaceColor',color7,'EdgeColor',color7); hold on
xticks([1 2 4 5 7 8 10 11 13 14 16 17]); xticklabels({'4','5','4','6','4','7','5','6','5','7','6','7'});
box off; title('kFpe'); ylim([0 0.2])

subplot(2,6,11)
bar(1,CP2.V1.R.B,'FaceColor',color4,'EdgeColor',color4); hold on
bar(2,CP2.V1.R.B,'FaceColor',color5,'EdgeColor',color5); hold on
bar(4,CP2.V2.R.B,'FaceColor',color4,'EdgeColor',color4); hold on
bar(5,CP2.V2.R.B,'FaceColor',color6,'EdgeColor',color6); hold on
bar(7,CP2.V3.R.B,'FaceColor',color4,'EdgeColor',color4); hold on
bar(8,CP2.V3.R.B,'FaceColor',color7,'EdgeColor',color7); hold on
bar(10,CP2.V4.R.B,'FaceColor',color5,'EdgeColor',color5); hold on
bar(11,CP2.V4.R.B,'FaceColor',color6,'EdgeColor',color6); hold on
bar(13,CP2.V5.R.B,'FaceColor',color5,'EdgeColor',color5); hold on
bar(14,CP2.V5.R.B,'FaceColor',color7,'EdgeColor',color7); hold on
bar(16,CP2.V6.R.B,'FaceColor',color6,'EdgeColor',color6); hold on
bar(17,CP2.V6.R.B,'FaceColor',color7,'EdgeColor',color7); hold on
xticks([1 2 4 5 7 8 10 11 13 14 16 17]); xticklabels({'4','5','4','6','4','7','5','6','5','7','6','7'});
box off; title('B'); ylim([0 0.02])

suptitle('CP2')
saveas(gcf,'C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/Figures/CP2_MT2.fig')
saveas(gcf,'C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/Figures/CP2_MT2.png')


%% CP 8
color2 = [0 48 73]./255;
color3 = [214 40 40]./255;
color4 = [247 127 0]./255;
color5 = [252 191 73]./255;

figure()
trials = [2 3; 2 4; 2 5;3 4; 3 5; 4 5]; 
for i = 1:length(trials)
    
    name  = ['T',num2str(trials(i,1)),'-T',num2str(trials(i,2))];
    Name  = ['CP8_', name,'_Opt8_MT.mat']
    fname = ['V',num2str(i)];
    CP8.(char(fname)) = load(Name);
    
    subplot(2,6,i)
    plot(CP8.(char(fname)).R.exp.qspline,'k','LineWidth',1.5); hold on
    plot(CP8.(char(fname)).R.x,'r','LineWidth',1.5); hold on
    title(name); box off;
    ylim([-2 0]); xlim([0 2000])
end

subplot(2,6,7)
bar(1,CP8.V1.R.a_ext(1),'FaceColor',color2,'EdgeColor',color2); hold on
bar(2,CP8.V1.R.a_ext(2),'FaceColor',color3,'EdgeColor',color3); hold on
bar(4,CP8.V2.R.a_ext(1),'FaceColor',color2,'EdgeColor',color2); hold on
bar(5,CP8.V2.R.a_ext(2),'FaceColor',color4,'EdgeColor',color4); hold on
bar(7,CP8.V3.R.a_ext(1),'FaceColor',color2,'EdgeColor',color2); hold on
bar(8,CP8.V3.R.a_ext(2),'FaceColor',color5,'EdgeColor',color5); hold on
bar(10,CP8.V4.R.a_ext(1),'FaceColor',color3,'EdgeColor',color3); hold on
bar(11,CP8.V4.R.a_ext(2),'FaceColor',color4,'EdgeColor',color4); hold on
bar(13,CP8.V5.R.a_ext(1),'FaceColor',color3,'EdgeColor',color3); hold on
bar(14,CP8.V5.R.a_ext(2),'FaceColor',color5,'EdgeColor',color5); hold on
bar(16,CP8.V6.R.a_ext(1),'FaceColor',color4,'EdgeColor',color4); hold on
bar(17,CP8.V6.R.a_ext(2),'FaceColor',color5,'EdgeColor',color5); hold on
xticks([1 2 4 5 7 8 10 11 13 14 16 17]); xticklabels({'2','3','2','4','2','5','3','4','3','5','4','5'});
box off; title('a ext'); ylim([0 0.02])

subplot(2,6,8)
bar(1,CP8.V1.R.a_flex(1),'FaceColor',color2,'EdgeColor',color2); hold on
bar(2,CP8.V1.R.a_flex(2),'FaceColor',color3,'EdgeColor',color3); hold on
bar(4,CP8.V2.R.a_flex(1),'FaceColor',color2,'EdgeColor',color2); hold on
bar(5,CP8.V2.R.a_flex(2),'FaceColor',color4,'EdgeColor',color4); hold on
bar(7,CP8.V3.R.a_flex(1),'FaceColor',color2,'EdgeColor',color2); hold on
bar(8,CP8.V3.R.a_flex(2),'FaceColor',color5,'EdgeColor',color5); hold on
bar(10,CP8.V4.R.a_flex(1),'FaceColor',color3,'EdgeColor',color3); hold on
bar(11,CP8.V4.R.a_flex(2),'FaceColor',color4,'EdgeColor',color4); hold on
bar(13,CP8.V5.R.a_flex(1),'FaceColor',color3,'EdgeColor',color3); hold on
bar(14,CP8.V5.R.a_flex(2),'FaceColor',color5,'EdgeColor',color5); hold on
bar(16,CP8.V6.R.a_flex(1),'FaceColor',color4,'EdgeColor',color4); hold on
bar(17,CP8.V6.R.a_flex(2),'FaceColor',color5,'EdgeColor',color5); hold on
xticks([1 2 4 5 7 8 10 11 13 14 16 17]); xticklabels({'2','3','2','4','2','5','3','4','3','5','4','5'});
box off; title('a flex'); ylim([0 0.02])

subplot(2,6,9)
bar(1,CP8.V1.R.kR(1),'FaceColor',color2,'EdgeColor',color2); hold on
bar(2,CP8.V1.R.kR(2),'FaceColor',color3,'EdgeColor',color3); hold on
bar(4,CP8.V2.R.kR(1),'FaceColor',color2,'EdgeColor',color2); hold on
bar(5,CP8.V2.R.kR(2),'FaceColor',color4,'EdgeColor',color4); hold on
bar(7,CP8.V3.R.kR(1),'FaceColor',color2,'EdgeColor',color2); hold on
bar(8,CP8.V3.R.kR(2),'FaceColor',color5,'EdgeColor',color5); hold on
bar(10,CP8.V4.R.kR(1),'FaceColor',color3,'EdgeColor',color3); hold on
bar(11,CP8.V4.R.kR(2),'FaceColor',color4,'EdgeColor',color4); hold on
bar(13,CP8.V5.R.kR(1),'FaceColor',color3,'EdgeColor',color3); hold on
bar(14,CP8.V5.R.kR(2),'FaceColor',color5,'EdgeColor',color5); hold on
bar(16,CP8.V6.R.kR(1),'FaceColor',color4,'EdgeColor',color4); hold on
bar(17,CP8.V6.R.kR(2),'FaceColor',color5,'EdgeColor',color5); hold on
xticks([1 2 4 5 7 8 10 11 13 14 16 17]); xticklabels({'2','3','2','4','2','5','3','4','3','5','4','5'});
box off; title('kR'); ylim([0 10])

subplot(2,6,10)
bar(1,CP8.V1.R.kFpe,'FaceColor',color2,'EdgeColor',color2); hold on
bar(2,CP8.V1.R.kFpe,'FaceColor',color3,'EdgeColor',color3); hold on
bar(4,CP8.V2.R.kFpe,'FaceColor',color2,'EdgeColor',color2); hold on
bar(5,CP8.V2.R.kFpe,'FaceColor',color4,'EdgeColor',color4); hold on
bar(7,CP8.V3.R.kFpe,'FaceColor',color2,'EdgeColor',color2); hold on
bar(8,CP8.V3.R.kFpe,'FaceColor',color5,'EdgeColor',color5); hold on
bar(10,CP8.V4.R.kFpe,'FaceColor',color3,'EdgeColor',color3); hold on
bar(11,CP8.V4.R.kFpe,'FaceColor',color4,'EdgeColor',color4); hold on
bar(13,CP8.V5.R.kFpe,'FaceColor',color3,'EdgeColor',color3); hold on
bar(14,CP8.V5.R.kFpe,'FaceColor',color5,'EdgeColor',color5); hold on
bar(16,CP8.V6.R.kFpe,'FaceColor',color4,'EdgeColor',color4); hold on
bar(17,CP8.V6.R.kFpe,'FaceColor',color5,'EdgeColor',color5); hold on
xticks([1 2 4 5 7 8 10 11 13 14 16 17]); xticklabels({'2','3','2','4','2','5','3','4','3','5','4','5'});
box off; title('kFpe');ylim([0 0.2])

subplot(2,6,11)
bar(1,CP8.V1.R.B,'FaceColor',color2,'EdgeColor',color2); hold on
bar(2,CP8.V1.R.B,'FaceColor',color3,'EdgeColor',color3); hold on
bar(4,CP8.V2.R.B,'FaceColor',color2,'EdgeColor',color2); hold on
bar(5,CP8.V2.R.B,'FaceColor',color4,'EdgeColor',color4); hold on
bar(7,CP8.V3.R.B,'FaceColor',color2,'EdgeColor',color2); hold on
bar(8,CP8.V3.R.B,'FaceColor',color5,'EdgeColor',color5); hold on
bar(10,CP8.V4.R.B,'FaceColor',color3,'EdgeColor',color3); hold on
bar(11,CP8.V4.R.B,'FaceColor',color4,'EdgeColor',color4); hold on
bar(13,CP8.V5.R.B,'FaceColor',color3,'EdgeColor',color3); hold on
bar(14,CP8.V5.R.B,'FaceColor',color5,'EdgeColor',color5); hold on
bar(16,CP8.V6.R.B,'FaceColor',color4,'EdgeColor',color4); hold on
bar(17,CP8.V6.R.B,'FaceColor',color5,'EdgeColor',color5); hold on
xticks([1 2 4 5 7 8 10 11 13 14 16 17]); xticklabels({'2','3','2','4','2','5','3','4','3','5','4','5'});
box off; title('B'); ylim([0 0.15])

suptitle('CP8')
saveas(gcf,'C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/Figures/CP8_MT2.fig')
saveas(gcf,'C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/Figures/CP8_MT2.png')


%% TD5
color1 = [0 48 73]./255;
color2 = [214 40 40]./255;

figure()
trials = [1 2]; 
for i = 1
    
    name  = ['T',num2str(trials(i,1)),'-T',num2str(trials(i,2))];
    Name  = ['TD5_', name,'_Opt8_MT.mat']
    fname = ['V',num2str(i)];
    TD5.(char(fname)) = load(Name);
    
    subplot(2,5,i)
    plot(TD5.(char(fname)).R.exp.qspline,'k','LineWidth',1.5); hold on
    plot(TD5.(char(fname)).R.x,'r','LineWidth',1.5); hold on
    title(name); box off;
    ylim([-2.5 0]); xlim([0 4000])
end

subplot(2,5,6)
bar(1,TD5.V1.R.a_ext(1),'FaceColor',color1,'EdgeColor',color1); hold on
bar(2,TD5.V1.R.a_ext(2),'FaceColor',color2,'EdgeColor',color2); hold on
xticks([1 2 ]); xticklabels({'1','2'});
box off; title('a ext'); ylim([0 0.02])

subplot(2,5,7)
bar(1,TD5.V1.R.a_flex(1),'FaceColor',color1,'EdgeColor',color1); hold on
bar(2,TD5.V1.R.a_flex(2),'FaceColor',color2,'EdgeColor',color2); hold on
xticks([1 2 ]); xticklabels({'1','2'});
box off; title('a flex'); ylim([0 0.02])

subplot(2,5,8)
bar(1,TD5.V1.R.kR(1),'FaceColor',color1,'EdgeColor',color1); hold on
bar(2,TD5.V1.R.kR(2),'FaceColor',color2,'EdgeColor',color2); hold on
xticks([1 2 ]); xticklabels({'1','2'});
box off; title('kR'); ylim([0 1])

subplot(2,5,9)
bar(1,TD5.V1.R.kFpe,'FaceColor',color1,'EdgeColor',color1); hold on
bar(2,TD5.V1.R.kFpe,'FaceColor',color2,'EdgeColor',color2); hold on
xticks([1 2 ]); xticklabels({'1','2'});
box off; title('kFpe');ylim([0 0.2])

subplot(2,5,10)
bar(1,TD5.V1.R.B,'FaceColor',color1,'EdgeColor',color1); hold on
bar(2,TD5.V1.R.B,'FaceColor',color2,'EdgeColor',color2); hold on
xticks([1 2 ]); xticklabels({'1','2'});
box off; title('B'); ylim([0 0.15])

suptitle('TD5')
saveas(gcf,'C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/Figures/TD5_MT2.fig')
saveas(gcf,'C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/Figures/TD5_MT2.png')


%% CP 4
color3 = [0 48 73]./255;
color4 = [214 40 40]./255;

figure()
trials = [3 4]; 
for i = 1
    
    name  = ['T',num2str(trials(i,1)),'-T',num2str(trials(i,2))];
    Name  = ['CP4_', name,'_Opt8_MT.mat']
    fname = ['V',num2str(i)];
    CP4.(char(fname)) = load(Name);
    
    subplot(2,5,i)
    plot(CP4.(char(fname)).R.exp.qspline,'k','LineWidth',1.5); hold on
    plot(CP4.(char(fname)).R.x,'r','LineWidth',1.5); hold on
    title(name); box off;
    ylim([-2.5 0]); xlim([0 4000])
end

subplot(2,5,6)
bar(1,CP4.V1.R.a_ext(1),'FaceColor',color3,'EdgeColor',color3); hold on
bar(2,CP4.V1.R.a_ext(2),'FaceColor',color4,'EdgeColor',color4); hold on
xticks([1 2 ]); xticklabels({'1','2'});
box off; title('a ext'); ylim([0 0.02])

subplot(2,5,7)
bar(1,CP4.V1.R.a_flex(1),'FaceColor',color3,'EdgeColor',color3); hold on
bar(2,CP4.V1.R.a_flex(2),'FaceColor',color4,'EdgeColor',color4); hold on
xticks([1 2 ]); xticklabels({'1','2'});
box off; title('a flex'); ylim([0 0.02])

subplot(2,5,8)
bar(1,CP4.V1.R.kR(1),'FaceColor',color3,'EdgeColor',color3); hold on
bar(2,CP4.V1.R.kR(2),'FaceColor',color4,'EdgeColor',color4); hold on
xticks([1 2 ]); xticklabels({'1','2'});
box off; title('kR'); ylim([0 1])

subplot(2,5,9)
bar(1,CP4.V1.R.kFpe,'FaceColor',color3,'EdgeColor',color3); hold on
bar(2,CP4.V1.R.kFpe,'FaceColor',color4,'EdgeColor',color4); hold on
xticks([1 2 ]); xticklabels({'1','2'});
box off; title('kFpe');ylim([0 0.2])

subplot(2,5,10)
bar(1,CP4.V1.R.B,'FaceColor',color3,'EdgeColor',color3); hold on
bar(2,CP4.V1.R.B,'FaceColor',color4,'EdgeColor',color4); hold on
xticks([1 2 ]); xticklabels({'1','2'});
box off; title('B'); ylim([0 0.15])

suptitle('CP4')
saveas(gcf,'C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/Figures/CP4_MT2.fig')
saveas(gcf,'C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/Figures/CP4_MT2.png')

%% CP 10
color12 = [0 48 73]./255;
color13 = [214 40 40]./255;
color14 = [247 127 0]./255;

figure()
trials = [12 14; 13 14]; 
for i = 1:length(trials)
    
    name  = ['T',num2str(trials(i,1)),'-T',num2str(trials(i,2))];
    Name  = ['CP10_', name,'_Opt8_MT.mat']
    fname = ['V',num2str(i)];
    CP10.(char(fname)) = load(Name);
    
    subplot(2,5,i)
    plot(CP10.(char(fname)).R.exp.qspline,'k','LineWidth',1.5); hold on
    plot(CP10.(char(fname)).R.x,'r','LineWidth',1.5); hold on
    title(name); box off;
    ylim([-2.5 0]); xlim([0 2000])
end

subplot(2,5,6)
bar(1,CP10.V1.R.a_ext(1),'FaceColor',color12,'EdgeColor',color12); hold on
bar(2,CP10.V1.R.a_ext(2),'FaceColor',color14,'EdgeColor',color14); hold on
bar(4,CP10.V2.R.a_ext(1),'FaceColor',color13,'EdgeColor',color13); hold on
bar(5,CP10.V2.R.a_ext(2),'FaceColor',color14,'EdgeColor',color14); hold on
xticks([1 2 4 5 ]); xticklabels({'12','13','13','14'});
box off; title('a ext'); ylim([0 0.03])

subplot(2,5,7)
bar(1,CP10.V1.R.a_flex(1),'FaceColor',color12,'EdgeColor',color12); hold on
bar(2,CP10.V1.R.a_flex(2),'FaceColor',color14,'EdgeColor',color14); hold on
bar(4,CP10.V2.R.a_flex(1),'FaceColor',color13,'EdgeColor',color13); hold on
bar(5,CP10.V2.R.a_flex(2),'FaceColor',color14,'EdgeColor',color14); hold on
xticks([1 2 4 5 ]); xticklabels({'12','13','13','14'});
box off; title('a flex'); ylim([0 0.03])

subplot(2,5,8)
bar(1,CP10.V1.R.kR(1),'FaceColor',color12,'EdgeColor',color12); hold on
bar(2,CP10.V1.R.kR(2),'FaceColor',color14,'EdgeColor',color14); hold on
bar(4,CP10.V2.R.kR(1),'FaceColor',color13,'EdgeColor',color13); hold on
bar(5,CP10.V2.R.kR(2),'FaceColor',color14,'EdgeColor',color14); hold on
xticks([1 2 4 5 ]); xticklabels({'12','13','13','14'});
box off; title('kR'); ylim([0 10])

subplot(2,5,9)
bar(1,CP10.V1.R.kFpe,'FaceColor',color12,'EdgeColor',color12); hold on
bar(2,CP10.V1.R.kFpe,'FaceColor',color14,'EdgeColor',color14); hold on
bar(4,CP10.V2.R.kFpe,'FaceColor',color13,'EdgeColor',color13); hold on
bar(5,CP10.V2.R.kFpe,'FaceColor',color14,'EdgeColor',color14); hold on
xticks([1 2 4 5 ]); xticklabels({'12','13','13','14'});
box off; title('kFpe');ylim([0 0.2])

subplot(2,5,10)
bar(1,CP10.V1.R.B,'FaceColor',color12,'EdgeColor',color12); hold on
bar(2,CP10.V1.R.B,'FaceColor',color14,'EdgeColor',color14); hold on
bar(4,CP10.V2.R.B,'FaceColor',color13,'EdgeColor',color13); hold on
bar(5,CP10.V2.R.B,'FaceColor',color14,'EdgeColor',color14); hold on
xticks([1 2 4 5 ]); xticklabels({'12','13','13','14'});
box off; title('B'); ylim([0 0.2])

suptitle('CP10')
saveas(gcf,'C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/Figures/CP10_MT2.fig')
saveas(gcf,'C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/Figures/CP10_MT2.png')

%% CP 11
color2 = [0 48 73]./255;
color3 = [214 40 40]./255;
color4 = [247 127 0]./255;

figure()
trials = [2 3; 2 4; 3 4;]; 
for i = 1:length(trials)
    
    name  = ['T',num2str(trials(i,1)),'-T',num2str(trials(i,2))];
    Name  = ['CP11_', name,'_Opt8_MT.mat']
    fname = ['V',num2str(i)];
    CP11.(char(fname)) = load(Name);
    
    subplot(2,5,i)
    plot(CP11.(char(fname)).R.exp.qspline,'k','LineWidth',1.5); hold on
    plot(CP11.(char(fname)).R.x,'r','LineWidth',1.5); hold on
    title(name); box off;
    ylim([-2 0]); xlim([0 2000])
end

subplot(2,5,6)
bar(1,CP11.V1.R.a_ext(1),'FaceColor',color2,'EdgeColor',color2); hold on
bar(2,CP11.V1.R.a_ext(2),'FaceColor',color3,'EdgeColor',color3); hold on
bar(4,CP11.V2.R.a_ext(1),'FaceColor',color2,'EdgeColor',color2); hold on
bar(5,CP11.V2.R.a_ext(2),'FaceColor',color4,'EdgeColor',color4); hold on
bar(7,CP11.V3.R.a_ext(1),'FaceColor',color3,'EdgeColor',color3); hold on
bar(8,CP11.V3.R.a_ext(2),'FaceColor',color4,'EdgeColor',color4); hold on
xticks([1 2 4 5 7 8 ]); xticklabels({'2','3','2','4','3','4'});
box off; title('a ext'); ylim([0 0.035])

subplot(2,5,7)
bar(1,CP11.V1.R.a_flex(1),'FaceColor',color2,'EdgeColor',color2); hold on
bar(2,CP11.V1.R.a_flex(2),'FaceColor',color3,'EdgeColor',color3); hold on
bar(4,CP11.V2.R.a_flex(1),'FaceColor',color2,'EdgeColor',color2); hold on
bar(5,CP11.V2.R.a_flex(2),'FaceColor',color4,'EdgeColor',color4); hold on
bar(7,CP11.V3.R.a_flex(1),'FaceColor',color3,'EdgeColor',color3); hold on
bar(8,CP11.V3.R.a_flex(2),'FaceColor',color4,'EdgeColor',color4); hold on
xticks([1 2 4 5 7 8 ]); xticklabels({'2','3','2','4','3','4'});
box off; title('a flex'); ylim([0 0.035])

subplot(2,5,8)
bar(1,CP11.V1.R.kR(1),'FaceColor',color2,'EdgeColor',color2); hold on
bar(2,CP11.V1.R.kR(2),'FaceColor',color3,'EdgeColor',color3); hold on
bar(4,CP11.V2.R.kR(1),'FaceColor',color2,'EdgeColor',color2); hold on
bar(5,CP11.V2.R.kR(2),'FaceColor',color4,'EdgeColor',color4); hold on
bar(7,CP11.V3.R.kR(1),'FaceColor',color3,'EdgeColor',color3); hold on
bar(8,CP11.V3.R.kR(2),'FaceColor',color4,'EdgeColor',color4); hold on
xticks([1 2 4 5 7 8 ]); xticklabels({'2','3','2','4','3','4'});
box off; title('kR'); ylim([0 10])

subplot(2,5,9)
bar(1,CP11.V1.R.kFpe,'FaceColor',color2,'EdgeColor',color2); hold on
bar(2,CP11.V1.R.kFpe,'FaceColor',color3,'EdgeColor',color3); hold on
bar(4,CP11.V2.R.kFpe,'FaceColor',color2,'EdgeColor',color2); hold on
bar(5,CP11.V2.R.kFpe,'FaceColor',color4,'EdgeColor',color4); hold on
bar(7,CP11.V3.R.kFpe,'FaceColor',color3,'EdgeColor',color3); hold on
bar(8,CP11.V3.R.kFpe,'FaceColor',color4,'EdgeColor',color4); hold on
xticks([1 2 4 5 7 8 ]); xticklabels({'2','3','2','4','3','4'});
box off; title('kFpe');ylim([0 0.2])

subplot(2,5,10)
bar(1,CP11.V1.R.B,'FaceColor',color2,'EdgeColor',color2); hold on
bar(2,CP11.V1.R.B,'FaceColor',color3,'EdgeColor',color3); hold on
bar(4,CP11.V2.R.B,'FaceColor',color2,'EdgeColor',color2); hold on
bar(5,CP11.V2.R.B,'FaceColor',color4,'EdgeColor',color4); hold on
bar(7,CP11.V3.R.B,'FaceColor',color3,'EdgeColor',color3); hold on
bar(8,CP11.V3.R.B,'FaceColor',color4,'EdgeColor',color4); hold on
xticks([1 2 4 5 7 8 ]); xticklabels({'2','3','2','4','3','4'});
box off; title('B'); ylim([0 0.15])

suptitle('CP11')
saveas(gcf,'C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/Figures/CP11_MT2.fig')
saveas(gcf,'C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/Figures/CP11_MT2.png')

%% CP 9
color2 = [0 48 73]./255;
color3 = [214 40 40]./255;
color4 = [247 127 0]./255;
color5 = [252 191 73]./255;
color6 = [234 226 183]./255;

figure()
trials = [2 3; 2 4; 2 5; 2 6; 3 4; 3 6; 4 5; 4 6; 5 6]; 
for i = 1:length(trials)
    
    name  = ['T',num2str(trials(i,1)),'-T',num2str(trials(i,2))];
    Name  = ['CP9_', name,'_Opt8_MT.mat']
    fname = ['V',num2str(i)];
    CP9.(char(fname)) = load(Name);
    
    subplot(2,10,i)
    plot(CP9.(char(fname)).R.exp.qspline,'k','LineWidth',1.5); hold on
    plot(CP9.(char(fname)).R.x,'r','LineWidth',1.5); hold on
    title(name); box off;
    ylim([-2.5 0]); xlim([0 3500])
end

subplot(2,9,10)
bar(1,CP9.V1.R.a_ext(1),'FaceColor',color2,'EdgeColor',color2); hold on
bar(2,CP9.V1.R.a_ext(2),'FaceColor',color3,'EdgeColor',color3); hold on
bar(4,CP9.V2.R.a_ext(1),'FaceColor',color2,'EdgeColor',color2); hold on
bar(5,CP9.V2.R.a_ext(2),'FaceColor',color4,'EdgeColor',color4); hold on
bar(7,CP9.V3.R.a_ext(1),'FaceColor',color2,'EdgeColor',color2); hold on
bar(8,CP9.V3.R.a_ext(2),'FaceColor',color5,'EdgeColor',color5); hold on
bar(10,CP9.V4.R.a_ext(1),'FaceColor',color2,'EdgeColor',color2); hold on
bar(11,CP9.V4.R.a_ext(2),'FaceColor',color6,'EdgeColor',color6); hold on
bar(13,CP9.V5.R.a_ext(1),'FaceColor',color3,'EdgeColor',color3); hold on
bar(14,CP9.V5.R.a_ext(2),'FaceColor',color4,'EdgeColor',color4); hold on
bar(16,CP9.V6.R.a_ext(1),'FaceColor',color3,'EdgeColor',color3); hold on
bar(17,CP9.V6.R.a_ext(2),'FaceColor',color6,'EdgeColor',color6); hold on
bar(19,CP9.V7.R.a_ext(1),'FaceColor',color4,'EdgeColor',color4); hold on
bar(20,CP9.V7.R.a_ext(2),'FaceColor',color5,'EdgeColor',color5); hold on
bar(22,CP9.V8.R.a_ext(1),'FaceColor',color4,'EdgeColor',color4); hold on
bar(23,CP9.V8.R.a_ext(2),'FaceColor',color6,'EdgeColor',color6); hold on
bar(25,CP9.V9.R.a_ext(1),'FaceColor',color5,'EdgeColor',color5); hold on
bar(26,CP9.V9.R.a_ext(2),'FaceColor',color6,'EdgeColor',color6); hold on
xticks([1 2 4 5 7 8 10 11 13 14 16 17 19 20 22 23 25 26]); xticklabels({'2','3','2','4','2','5','2','6','3','4','3','6','4','5','4','6','5','6'});
box off; title('a ext'); ylim([0 0.01])

subplot(2,9,11)
bar(1,CP9.V1.R.a_flex(1),'FaceColor',color2,'EdgeColor',color2); hold on
bar(2,CP9.V1.R.a_flex(2),'FaceColor',color3,'EdgeColor',color3); hold on
bar(4,CP9.V2.R.a_flex(1),'FaceColor',color2,'EdgeColor',color2); hold on
bar(5,CP9.V2.R.a_flex(2),'FaceColor',color4,'EdgeColor',color4); hold on
bar(7,CP9.V3.R.a_flex(1),'FaceColor',color2,'EdgeColor',color2); hold on
bar(8,CP9.V3.R.a_flex(2),'FaceColor',color5,'EdgeColor',color5); hold on
bar(10,CP9.V4.R.a_flex(1),'FaceColor',color2,'EdgeColor',color2); hold on
bar(11,CP9.V4.R.a_flex(2),'FaceColor',color6,'EdgeColor',color6); hold on
bar(13,CP9.V5.R.a_flex(1),'FaceColor',color3,'EdgeColor',color3); hold on
bar(14,CP9.V5.R.a_flex(2),'FaceColor',color4,'EdgeColor',color4); hold on
bar(16,CP9.V6.R.a_flex(1),'FaceColor',color3,'EdgeColor',color3); hold on
bar(17,CP9.V6.R.a_flex(2),'FaceColor',color6,'EdgeColor',color6); hold on
bar(19,CP9.V7.R.a_flex(1),'FaceColor',color4,'EdgeColor',color4); hold on
bar(20,CP9.V7.R.a_flex(2),'FaceColor',color5,'EdgeColor',color5); hold on
bar(22,CP9.V8.R.a_flex(1),'FaceColor',color4,'EdgeColor',color4); hold on
bar(23,CP9.V8.R.a_flex(2),'FaceColor',color6,'EdgeColor',color6); hold on
bar(25,CP9.V9.R.a_flex(1),'FaceColor',color5,'EdgeColor',color5); hold on
bar(26,CP9.V9.R.a_flex(2),'FaceColor',color6,'EdgeColor',color6); hold on
xticks([1 2 4 5 7 8 10 11 13 14 16 17 19 20 22 23 25 26]); xticklabels({'2','3','2','4','2','5','2','6','3','4','3','6','4','5','4','6','5','6'});
box off; title('a flex'); ylim([0 0.01])

subplot(2,9,12)
bar(1,CP9.V1.R.kR(1),'FaceColor',color2,'EdgeColor',color2); hold on
bar(2,CP9.V1.R.kR(2),'FaceColor',color3,'EdgeColor',color3); hold on
bar(4,CP9.V2.R.kR(1),'FaceColor',color2,'EdgeColor',color2); hold on
bar(5,CP9.V2.R.kR(2),'FaceColor',color4,'EdgeColor',color4); hold on
bar(7,CP9.V3.R.kR(1),'FaceColor',color2,'EdgeColor',color2); hold on
bar(8,CP9.V3.R.kR(2),'FaceColor',color5,'EdgeColor',color5); hold on
bar(10,CP9.V4.R.kR(1),'FaceColor',color2,'EdgeColor',color2); hold on
bar(11,CP9.V4.R.kR(2),'FaceColor',color6,'EdgeColor',color6); hold on
bar(13,CP9.V5.R.kR(1),'FaceColor',color3,'EdgeColor',color3); hold on
bar(14,CP9.V5.R.kR(2),'FaceColor',color4,'EdgeColor',color4); hold on
bar(16,CP9.V6.R.kR(1),'FaceColor',color3,'EdgeColor',color3); hold on
bar(17,CP9.V6.R.kR(2),'FaceColor',color6,'EdgeColor',color6); hold on
bar(19,CP9.V7.R.kR(1),'FaceColor',color4,'EdgeColor',color4); hold on
bar(20,CP9.V7.R.kR(2),'FaceColor',color5,'EdgeColor',color5); hold on
bar(22,CP9.V8.R.kR(1),'FaceColor',color4,'EdgeColor',color4); hold on
bar(23,CP9.V8.R.kR(2),'FaceColor',color6,'EdgeColor',color6); hold on
bar(25,CP9.V9.R.kR(1),'FaceColor',color5,'EdgeColor',color5); hold on
bar(26,CP9.V9.R.kR(2),'FaceColor',color6,'EdgeColor',color6); hold on
xticks([1 2 4 5 7 8 10 11 13 14 16 17 19 20 22 23 25 26]); xticklabels({'2','3','2','4','2','5','2','6','3','4','3','6','4','5','4','6','5','6'});
box off; title('kR'); ylim([0 1])

subplot(2,9,13)
bar(1,CP9.V1.R.kFpe,'FaceColor',color2,'EdgeColor',color2); hold on
bar(2,CP9.V1.R.kFpe,'FaceColor',color3,'EdgeColor',color3); hold on
bar(4,CP9.V2.R.kFpe,'FaceColor',color2,'EdgeColor',color2); hold on
bar(5,CP9.V2.R.kFpe,'FaceColor',color4,'EdgeColor',color4); hold on
bar(7,CP9.V3.R.kFpe,'FaceColor',color2,'EdgeColor',color2); hold on
bar(8,CP9.V3.R.kFpe,'FaceColor',color5,'EdgeColor',color5); hold on
bar(10,CP9.V4.R.kFpe,'FaceColor',color2,'EdgeColor',color2); hold on
bar(11,CP9.V4.R.kFpe,'FaceColor',color6,'EdgeColor',color6); hold on
bar(13,CP9.V5.R.kFpe,'FaceColor',color3,'EdgeColor',color3); hold on
bar(14,CP9.V5.R.kFpe,'FaceColor',color4,'EdgeColor',color4); hold on
bar(16,CP9.V6.R.kFpe,'FaceColor',color3,'EdgeColor',color3); hold on
bar(17,CP9.V6.R.kFpe,'FaceColor',color6,'EdgeColor',color6); hold on
bar(19,CP9.V7.R.kFpe,'FaceColor',color4,'EdgeColor',color4); hold on
bar(20,CP9.V7.R.kFpe,'FaceColor',color5,'EdgeColor',color5); hold on
bar(22,CP9.V8.R.kFpe,'FaceColor',color4,'EdgeColor',color4); hold on
bar(23,CP9.V8.R.kFpe,'FaceColor',color6,'EdgeColor',color6); hold on
bar(25,CP9.V9.R.kFpe,'FaceColor',color5,'EdgeColor',color5); hold on
bar(26,CP9.V9.R.kFpe,'FaceColor',color6,'EdgeColor',color6); hold on
xticks([1 2 4 5 7 8 10 11 13 14 16 17 19 20 22 23 25 26]); xticklabels({'2','3','2','4','2','5','2','6','3','4','3','6','4','5','4','6','5','6'});
box off; title('kFpe');ylim([0 0.2])

subplot(2,9,14)
bar(1,CP9.V1.R.B,'FaceColor',color2,'EdgeColor',color2); hold on
bar(2,CP9.V1.R.B,'FaceColor',color3,'EdgeColor',color3); hold on
bar(4,CP9.V2.R.B,'FaceColor',color2,'EdgeColor',color2); hold on
bar(5,CP9.V2.R.B,'FaceColor',color4,'EdgeColor',color4); hold on
bar(7,CP9.V3.R.B,'FaceColor',color2,'EdgeColor',color2); hold on
bar(8,CP9.V3.R.B,'FaceColor',color5,'EdgeColor',color5); hold on
bar(10,CP9.V4.R.B,'FaceColor',color2,'EdgeColor',color2); hold on
bar(11,CP9.V4.R.B,'FaceColor',color6,'EdgeColor',color6); hold on
bar(13,CP9.V5.R.B,'FaceColor',color3,'EdgeColor',color3); hold on
bar(14,CP9.V5.R.B,'FaceColor',color4,'EdgeColor',color4); hold on
bar(16,CP9.V6.R.B,'FaceColor',color3,'EdgeColor',color3); hold on
bar(17,CP9.V6.R.B,'FaceColor',color6,'EdgeColor',color6); hold on
bar(19,CP9.V7.R.B,'FaceColor',color4,'EdgeColor',color4); hold on
bar(20,CP9.V7.R.B,'FaceColor',color5,'EdgeColor',color5); hold on
bar(22,CP9.V8.R.B,'FaceColor',color4,'EdgeColor',color4); hold on
bar(23,CP9.V8.R.B,'FaceColor',color6,'EdgeColor',color6); hold on
bar(25,CP9.V9.R.B,'FaceColor',color5,'EdgeColor',color5); hold on
bar(26,CP9.V9.R.B,'FaceColor',color6,'EdgeColor',color6); hold on
xticks([1 2 4 5 7 8 10 11 13 14 16 17 19 20 22 23 25 26]); xticklabels({'2','3','2','4','2','5','2','6','3','4','3','6','4','5','4','6','5','6'});
box off; title('B'); ylim([0 0.05])

suptitle('CP9')
saveas(gcf,'C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/Figures/CP9_MT2.fig')
saveas(gcf,'C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/Figures/CP9_MT2.png')

