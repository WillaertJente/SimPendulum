%% Plot to compare different versions 
% Jente Willaert - 14032022
color1 = [0 48 73]./255;
color2 = [69 123 157]./255;
color1 = [250 163 7]./255;

v1 = load('C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/SDR2_Post_T1_Opt7.mat'); 
v2 = load('C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/SDR2_Post_T1_Opt9_Yank.mat'); 

figure()
subplot(2,2,1)
plot(v1.R.exp.qspline,'k','LineWidth',1.5); hold on
plot(v1.R.x,'r','LineWidth',1.5); hold on
box off; title('Reflex on SRS')
subplot(2,2,2)
plot(v2.R.exp.qspline,'k','LineWidth',1.5); hold on
plot(v2.R.x,'r','LineWidth',1.5); hold on
box off; title('Reflex on force and Yank')

subplot(2,6,7,'replace')
bar(1,v1.R.a_ext,'FaceColor',color1,'EdgeColor',color1); hold on
bar(2,v2.R.a_ext,'FaceColor',color2,'EdgeColor',color2); hold on
title('a Ext'); box off; xticks([1 2]); xticklabels({'v1','v2'})

subplot(2,6,8)
bar(1,v1.R.a_flex,'FaceColor',color1,'EdgeColor',color1); hold on
bar(2,v2.R.a_flex,'FaceColor',color2,'EdgeColor',color2); hold on
title('a Flex'); box off; xticks([1 2]); xticklabels({'v1','v2'})

subplot(2,6,9)
bar(1,v1.R.kR,'FaceColor',color1,'EdgeColor',color1); hold on
bar(2,v2.R.kR,'FaceColor',color2,'EdgeColor',color2); hold on
title('R'); box off; xticks([1 2]); xticklabels({'v1','v2'})

subplot(2,6,10)
bar(1,v1.R.kFpe,'FaceColor',color1,'EdgeColor',color1); hold on
bar(2,v2.R.kFpe,'FaceColor',color2,'EdgeColor',color2); hold on
title('kFpe'); box off;xticks([1 2]); xticklabels({'v1','v2'})

subplot(2,6,11)
bar(1,v1.R.B,'FaceColor',color1,'EdgeColor',color1); hold on
bar(2,v2.R.B,'FaceColor',color2,'EdgeColor',color2); hold on
title('B'); box off; xticks([1 2]); xticklabels({'v1','v2'})

subplot(2,6,12)
bar(1,v1.R.J,'FaceColor',color1,'EdgeColor',color1); hold on
bar(2,v2.R.J,'FaceColor',color2,'EdgeColor',color2); hold on
box off; title('J'); xticks([1 2]); xticklabels({'v1','v2'})

suptitle('SDR2_Post_T1')