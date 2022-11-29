%Plot to compare different versions 
%Jente Willaert - 14032022

clear all; close all; clc 

color1 = [0 48 73]./255;
color2 = [69 123 157]./255;
color1 = [250 163 7]./255;

v1 = load(['C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/CP4_T3_OpT14_Activation_TanH_tresh_FS_IG3.mat']); 
v2 = load('C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/CP4_T3_OpT14_Activation_TanH_tresh_FS_IG3_testInertie2.mat'); 
v3 = load('C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/CP4_T3_OpT14_Activation_TanH_tresh_FS_IG3_testInertie3.mat'); 
%v4 = load('C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/CP4_T2_OpT24_Activation_TanH_tresh_FS2_IG4.mat'); 
%v5 = load('C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/CP4_T2_OpT24_Activation_TanH_tresh_FS2_IG5.mat'); 
%v6 = load('C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/CP4_T2_OpT24_Activation_TanH_tresh_FS2_IG6.mat'); 
%v7 = load('C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/CP4_T2_OpT24_Activation_TanH_tresh_FS2_IG7.mat'); 
%v8 = load('C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/CP4_T2_OpT24_Activation_TanH_tresh_FS2_IG8.mat'); 
%v9 = load('C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/CP4_T2_OpT24_Activation_TanH_tresh_FS2_IG9.mat'); 
%v10 = load('C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/CP4_T2_OpT24_Activation_TanH_tresh_FS2_IG10.mat'); 

figure()
subplot(2,2,1)
plot(v1.R.exp.qspline,'k','LineWidth',1.5); hold on
plot(v1.R.x,'Color',[0, 0.4470, 0.7410],'LineWidth',1.5); hold on
plot(v2.R.x,'Color',[0.8500, 0.3250, 0.0980],'LineWidth',1.5); hold on
plot(v3.R.x,'Color',[0.9290, 0.6940, 0.1250],'LineWidth',1.5); hold on
%plot(v4.R.x,'Color',[0.4940, 0.1840, 0.5560],'LineWidth',1.5); hold on
%plot(v5.R.x,'Color',[0.4660, 0.6740, 0.1880],'LineWidth',1.5); hold on
box off; title('O14 - Inertie 4 - inertie 5')
%subplot(2,2,2)
%plot(v2.R.exp.qspline,'k','LineWidth',1.5); hold on
%plot(v2.R.x,'Color',[0.8500, 0.3250, 0.0980],'LineWidth',1.5); hold on
%plot(v8.R.exp.qspline,'k','LineWidth',1.5); hold on
%plot(v6.R.x,'Color',[0.3010, 0.7450, 0.9330],'LineWidth',1.5); hold on
%plot(v7.R.x,'Color',[0.6350, 0.0780, 0.1840],'LineWidth',1.5); hold on
%plot(v8.R.x,'Color',[0.8, 0.4, 0],'LineWidth',1.5); hold on
%plot(v9.R.x,'Color',[0.9608    0.7922    0.7647],'LineWidth',1.5); hold on
%plot(v10.R.x,'Color',[0.5961    0.7569    0.8510],'LineWidth',1.5); hold on
box off; title('FS 100')

subplot(2,6,7,'replace')
bar(1,v1.R.a_ext,'FaceColor',[0, 0.4470, 0.7410]); hold on %,'FaceColor',color1,'EdgeColor',color1); hold on
bar(2,v2.R.a_ext,'FaceColor',[0.8500, 0.3250, 0.0980]); hold on %,'FaceColor',color2,'EdgeColor',color2); hold on
bar(3,v3.R.a_ext,'FaceColor',[0.9290, 0.6940, 0.1250]); hold on %,'FaceColor',color1,'EdgeColor',color1); hold on
%bar(4,v4.R.a_ext,'FaceColor',[0.4940, 0.1840, 0.5560]); hold on %,'FaceColor',color2,'EdgeColor',color2); hold on
%bar(5,v5.R.a_ext,'FaceColor',[0.4660, 0.6740, 0.1880]); hold on
%bar(6,v6.R.a_ext,'FaceColor',[0.3010, 0.7450, 0.9330]); hold on
%bar(7,v7.R.a_ext,'FaceColor',[0.6350, 0.0780, 0.1840]); hold on
%bar(8,v8.R.a_ext,'FaceColor',[0.8, 0.4, 0]); hold on
%bar(9,v9.R.a_ext,'FaceColor',[0.9608    0.7922    0.7647]); hold on
%bar(10,v10.R.a_ext,'FaceColor',[0.5961    0.7569    0.8510]); hold on
title('a Ext'); box off; xticks([1 2 3 4 5 6 7 8 9 10]); 

subplot(2,6,8)
bar(1,v1.R.a_flex,'FaceColor',[0, 0.4470, 0.7410]'); hold on %,'FaceColor',color1,'EdgeColor',color1); hold on
bar(2,v2.R.a_flex,'FaceColor',[0.8500, 0.3250, 0.0980]); hold on %,'FaceColor',color2,'EdgeColor',color2); hold on
bar(3,v3.R.a_flex,'FaceColor',[0.9290, 0.6940, 0.1250]); hold on
%bar(4,v4.R.a_flex,'FaceColor',[0.4940, 0.1840, 0.5560]); hold on
%bar(5,v5.R.a_flex,'FaceColor',[0.4660, 0.6740, 0.1880]); hold on
%bar(6,v6.R.a_flex,'FaceColor',[0.3010, 0.7450, 0.9330]); hold on
%bar(7,v7.R.a_flex,'FaceColor',[0.6350, 0.0780, 0.1840]); hold on
%bar(8,v8.R.a_flex,'FaceColor',[0.8, 0.4, 0]); hold on
%bar(9,v9.R.a_flex,'FaceColor',[0.9608    0.7922    0.7647]); hold on
%bar(10,v10.R.a_flex,'FaceColor',[0.5961    0.7569    0.8510]); hold on
title('a Flex'); box off; xticks([1 2 3 4 5 6 7 8 9 10]); 

subplot(2,6,9)
bar(1,v1.R.kR,'FaceColor',[0, 0.4470, 0.7410]); hold on %,'FaceColor',color1,'EdgeColor',color1); hold on
bar(2,v2.R.kR,'FaceColor',[0.8500, 0.3250, 0.0980]); hold on %,'FaceColor',color2,'EdgeColor',color2); hold on
bar(3,v3.R.kR,'FaceColor',[0.9290, 0.6940, 0.1250]); hold on %
%bar(4,v4.R.kR,'FaceColor',[0.4940, 0.1840, 0.5560]); hold on %
%bar(5,v5.R.kR,'FaceColor',[0.4660, 0.6740, 0.1880]); hold on %
%bar(6,v6.R.kR,'FaceColor',[0.3010, 0.7450, 0.9330]); hold on %
%bar(7,v7.R.kR,'FaceColor',[0.6350, 0.0780, 0.1840]); hold on %
%bar(8,v8.R.kR,'FaceColor',[0.8, 0.4, 0]); hold on %
%bar(9,v9.R.kR,'FaceColor',[0.9608    0.7922    0.7647]); hold on %
%bar(10,v10.R.kR,'FaceColor',[0.5961    0.7569    0.8510]); hold on %
title('R'); box off; xticks([1 2 3 4 5 6 7 8 9 10]); 

subplot(2,6,10)
bar(1,v1.R.kFpe,'FaceColor',[0, 0.4470, 0.7410]); hold on %,'FaceColor',color1,'EdgeColor',color1); hold on
bar(2,v2.R.kFpe,'FaceColor',[0.8500, 0.3250, 0.0980]); hold on %,'FaceColor',color2,'EdgeColor',color2); hold on
bar(3,v3.R.kFpe,'FaceColor',[0.9290, 0.6940, 0.1250]); hold on
%bar(4,v4.R.kFpe,'FaceColor',[0.4940, 0.1840, 0.5560]); hold on
%bar(5,v5.R.kFpe,'FaceColor',[0.4660, 0.6740, 0.1880]); hold on
%bar(6,v6.R.kFpe,'FaceColor',[0.3010, 0.7450, 0.9330]); hold on
%bar(7,v7.R.kFpe,'FaceColor',[0.6350, 0.0780, 0.1840]); hold on
%bar(8,v8.R.kFpe,'FaceColor',[0.8, 0.4, 0]); hold on
%bar(9,v9.R.kFpe,'FaceColor',[0.9608    0.7922    0.7647]); hold on
%bar(10,v10.R.kFpe,'FaceColor',[0.5961    0.7569    0.8510]); hold on
title('kFpe'); box off; xticks([1 2 3 4 5 6 7 8 9 10]);  

subplot(2,6,11)
bar(1,v1.R.B,'FaceColor',[0, 0.4470, 0.7410]); hold on %,'FaceColor',color1,'EdgeColor',color1); hold on
bar(2,v2.R.B,'FaceColor',[0.8500, 0.3250, 0.0980]); hold on %,'FaceColor',color2,'EdgeColor',color2); hold on
bar(3,v3.R.B,'FaceColor',[0.9290, 0.6940, 0.1250]); hold on
%bar(4,v4.R.B,'FaceColor',[0.4940, 0.1840, 0.5560]); hold on
%bar(5,v5.R.B,'FaceColor',[0.4660, 0.6740, 0.1880]); hold on
%bar(6,v6.R.B,'FaceColor',[0.3010, 0.7450, 0.9330]); hold on
%bar(7,v7.R.B,'FaceColor',[0.6350, 0.0780, 0.1840]); hold on
%bar(8,v8.R.B,'FaceColor',[0.8, 0.4, 0]); hold on
%bar(9,v9.R.B,'FaceColor',[0.9608    0.7922    0.7647]); hold on
%bar(10,v10.R.B,'FaceColor',[0.5961    0.7569    0.8510]); hold on
title('B'); box off; xticks([1 2 3 4 5 6 7 8 9 10]); 

subplot(2,6,12)
bar(1,v1.R.J,'FaceColor',[0, 0.4470, 0.7410]); hold on %,'FaceColor',color1,'EdgeColor',color1); hold on
bar(2,v2.R.J,'FaceColor',[0.8500, 0.3250, 0.0980]	); hold on %,'FaceColor',color2,'EdgeColor',color2); hold on
bar(3,v3.R.J,'FaceColor',[0.9290, 0.6940, 0.1250]); hold on
%bar(4,v4.R.J,'FaceColor',[0.4940, 0.1840, 0.5560]); hold on
%bar(5,v5.R.J,'FaceColor',[0.4660, 0.6740, 0.1880]); hold on
%bar(6,v6.R.J,'FaceColor',[0.3010, 0.7450, 0.9330]); hold on
%bar(7,v7.R.J,'FaceColor',[0.6350, 0.0780, 0.1840]); hold on
%bar(8,v8.R.J,'FaceColor',[0.8, 0.4, 0]); hold on
%bar(9,v9.R.J,'FaceColor',[0.9608    0.7922    0.7647]); hold on
%bar(10,v10.R.J,'FaceColor',[0.5961    0.7569    0.8510]); hold on
box off; title('J'); xticks([1 2 3 4 5 6 7 8 9 10]); 

sgtitle('CP4 - T3')