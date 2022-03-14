%% Plot to compare different versions 
% Jente Willaert - 14032022
color1 = [0 48 73]./255;
color2 = [214 40 40]./255;
color3 = [247 127 0]./255;

v1 = load('C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/CP2_T5_Opt7_ingekort_Bopt.mat'); 
v2 = load('C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/CP2_T5_Opt7_Bboundslower.mat'); 
v3 = load('C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/CP2_T5_Opt7_InitValue.mat'); 

figure()
subplot(2,6,1)
plot(v1.R.exp.qspline,'k','LineWidth',1.5); hold on
plot(v1.R.x,'r','LineWidth',1.5); hold on
box off; title('B Lb 0.01')
subplot(2,6,3)
plot(v2.R.exp.qspline,'k','LineWidth',1.5); hold on
plot(v2.R.x,'r','LineWidth',1.5); hold on
box off; title('B Lb 0.001 - local optima')
subplot(2,6,5)
plot(v3.R.exp.qspline,'k','LineWidth',1.5); hold on
plot(v3.R.x,'r','LineWidth',1.5); hold on
box off; title('B Lb 0.001 - init value v1')

subplot(2,6,7)
bar(1,v1.R.a_ext,'FaceColor',color1,'EdgeColor',color1); hold on
bar(2,v2.R.a_ext,'FaceColor',color2,'EdgeColor',color2); hold on
bar(3,v3.R.a_ext,'FaceColor',color3,'EdgeColor',color3); hold on
title('a Ext'); box off; 

subplot(2,6,8)
bar(1,v1.R.a_flex,'FaceColor',color1,'EdgeColor',color1); hold on
bar(2,v2.R.a_flex,'FaceColor',color2,'EdgeColor',color2); hold on
bar(3,v3.R.a_flex,'FaceColor',color3,'EdgeColor',color3); hold on
title('a Flex'); box off; 

subplot(2,6,9)
bar(1,v1.R.kR,'FaceColor',color1,'EdgeColor',color1); hold on
bar(2,v2.R.kR,'FaceColor',color2,'EdgeColor',color2); hold on
bar(3,v3.R.kR,'FaceColor',color3,'EdgeColor',color3); hold on
title('R'); box off; 

subplot(2,6,10)
bar(1,v1.R.kFpe,'FaceColor',color1,'EdgeColor',color1); hold on
bar(2,v2.R.kFpe,'FaceColor',color2,'EdgeColor',color2); hold on
bar(3,v3.R.kFpe,'FaceColor',color3,'EdgeColor',color3); hold on
title('kFpe'); box off;

subplot(2,6,11)
bar(1,v1.R.B,'FaceColor',color1,'EdgeColor',color1); hold on
bar(2,v2.R.B,'FaceColor',color2,'EdgeColor',color2); hold on
bar(3,v3.R.B,'FaceColor',color3,'EdgeColor',color3); hold on
title('B'); box off; 

subplot(2,6,12)
bar(1,v1.R.J,'FaceColor',color1,'EdgeColor',color1); hold on
bar(2,v2.R.J,'FaceColor',color2,'EdgeColor',color2); hold on
bar(3,v3.R.J,'FaceColor',color3,'EdgeColor',color3); hold on
box off; title('J')