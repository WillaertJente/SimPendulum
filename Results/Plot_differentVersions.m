% Plot to compare different versions 
%Jente Willaert - 14032022

clear all; close all; clc 

% color1 = [0 48 73]./255;
color2 = [69 123 157]./255;
color1 = [250 163 7]./255;

v1 = load('C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/CP9_T6_Opt12_IG_set4_BOunds_Set1_kost.mat'); 
v2 = load('C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/CP9_T6_Opt12_IG_set4_Bounds_Set2_kost.mat'); 
v3 = load('C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/CP9_T6_Opt12_IG_set4_Bounds_Set3_kost.mat'); 
v4 = load('C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/CP9_T6_Opt12_IG_set8_Bounds_Set1_kost.mat'); 
v5 = load('C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/CP9_T6_Opt12_IG_set8_Bounds_Set2_kost.mat'); 
v6 = load('C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/CP9_T6_Opt12_IG_set8_Bounds_Set3_kost.mat'); 
%v7 = load('C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/CP9_T4_Opt12_IG_set7.mat'); 
%v8 = load('C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/CP9_T1_Opt12_IG_set8.mat'); 
%v9 = load('C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/CP9_T1_Opt12_IG_set9.mat'); 
%v10 = load('C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/CP9_T1_Opt12_IG_set10.mat'); 

figure()
subplot(2,2,1)
plot(v1.R.exp.qspline,'k','LineWidth',1.5); hold on
plot(v1.R.x,'LineWidth',1.5); hold on
plot(v2.R.x,'LineWidth',1.5); hold on
plot(v3.R.x,'LineWidth',1.5); hold on
plot(v4.R.x,'LineWidth',1.5); hold on
plot(v5.R.x,'LineWidth',1.5); hold on
% box off; title('Same IG - Different ')
% subplot(2,2,2)
%plot(v4.R.exp.qspline,'k','LineWidth',1.5); hold on
% plot(v4.R.x,'LineWidth',1.5); hold on
% plot(v5.R.x,'LineWidth',1.5); hold on
plot(v6.R.x,'LineWidth',1.5); hold on
%plot(v9.R.x,'LineWidth',1.5); hold on
%plot(v10.R.x,'LineWidth',1.5); hold on
box off; title('Different bounds')

subplot(2,6,7,'replace')
bar(1,v1.R.a_ext); hold on %,'FaceColor',color1,'EdgeColor',color1); hold on
bar(2,v2.R.a_ext); hold on %,'FaceColor',color2,'EdgeColor',color2); hold on
bar(3,v3.R.a_ext); hold on %,'FaceColor',color1,'EdgeColor',color1); hold on
bar(4,v4.R.a_ext); hold on %,'FaceColor',color2,'EdgeColor',color2); hold on
bar(5,v5.R.a_ext); hold on
bar(6,v6.R.a_ext); hold on
%bar(7,v7.R.a_ext); hold on
%bar(8,v8.R.a_ext); hold on
%bar(9,v9.R.a_ext); hold on
%bar(10,v10.R.a_ext); hold on
title('a Ext'); box off; xticks([1 2 3 4 5 6])% 7 8 9 10]); 

subplot(2,6,8)
bar(1,v1.R.a_flex); hold on %,'FaceColor',color1,'EdgeColor',color1); hold on
bar(2,v2.R.a_flex); hold on %,'FaceColor',color2,'EdgeColor',color2); hold on
bar(3,v3.R.a_flex); hold on
bar(4,v4.R.a_flex); hold on
bar(5,v5.R.a_flex); hold on
bar(6,v6.R.a_flex); hold on
%bar(7,v7.R.a_flex); hold on
%bar(8,v8.R.a_flex); hold on
%bar(9,v9.R.a_flex); hold on
%bar(10,v10.R.a_flex); hold on
title('a Flex'); box off; xticks([1 2 3 4 5 6])% 7 8 9 10]); 

subplot(2,6,9)
bar(1,v1.R.kR); hold on %,'FaceColor',color1,'EdgeColor',color1); hold on
bar(2,v2.R.kR); hold on %,'FaceColor',color2,'EdgeColor',color2); hold on
bar(3,v3.R.kR); hold on %
bar(4,v4.R.kR); hold on %
bar(5,v5.R.kR); hold on %
bar(6,v6.R.kR); hold on %
%bar(7,v7.R.kR); hold on %
%bar(8,v8.R.kR); hold on %
%bar(9,v9.R.kR); hold on %
%bar(10,v10.R.kR); hold on %
title('R'); box off; xticks([1 2 3 4 5 6 ])%7 8 9 10]); 

subplot(2,6,10)
bar(1,v1.R.kFpe); hold on %,'FaceColor',color1,'EdgeColor',color1); hold on
bar(2,v2.R.kFpe); hold on %,'FaceColor',color2,'EdgeColor',color2); hold on
bar(3,v3.R.kFpe); hold on
%bar(4,v4.R.kFpe); hold on
%bar(5,v5.R.kFpe); hold on
%bar(6,v6.R.kFpe); hold on
%bar(7,v7.R.kFpe); hold on
%bar(8,v8.R.kFpe); hold on
%bar(9,v9.R.kFpe); hold on
%bar(10,v10.R.kFpe); hold on
title('kFpe'); box off; xticks([1 2 3 4 5 6])% 7 8 9 10]); 

subplot(2,6,11)
bar(1,v1.R.B); hold on %,'FaceColor',color1,'EdgeColor',color1); hold on
bar(2,v2.R.B); hold on %,'FaceColor',color2,'EdgeColor',color2); hold on
bar(3,v3.R.B); hold on
bar(4,v4.R.B); hold on
bar(5,v5.R.B); hold on
bar(6,v6.R.B); hold on
%bar(7,v7.R.B); hold on
%bar(8,v8.R.B); hold on
%bar(9,v9.R.B); hold on
%bar(10,v10.R.B); hold on
title('B'); box off; xticks([1 2 3 4 5 6])% 7 8 9 10]); 

subplot(2,6,12)
bar(1,v1.R.J); hold on %,'FaceColor',color1,'EdgeColor',color1); hold on
bar(2,v2.R.J); hold on %,'FaceColor',color2,'EdgeColor',color2); hold on
bar(3,v3.R.J); hold on
bar(4,v4.R.J); hold on
bar(5,v5.R.J); hold on
bar(6,v6.R.J); hold on
%bar(7,v7.R.J); hold on
%bar(8,v8.R.J); hold on
%bar(9,v9.R.J); hold on
%bar(10,v10.R.J); hold on
box off; title('J'); xticks([1 2 3 4 5 6 ])%7 8 9 10]); 

suptitle('CP9 - T6')