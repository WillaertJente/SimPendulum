%% Plot Pareto 

v1 = load('C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/CP11_T3_Opt7_wkR_0.1.mat');
v2 = load('C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/CP11_T3_Opt7_wkR_0.01.mat');
v3 = load('C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/CP11_T3_Opt7_wkR_0.001.mat');
v4 = load('C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/CP11_T3_Opt7_wkR_0.2.mat');
v5 = load('C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/CP11_T3_Opt7_wkR_0.02.mat');
v6 = load('C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/CP11_T3_Opt7_wkR_0.002.mat');
v7 = load('C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/CP11_T3_Opt7_wkR_0.5.mat');
v8 = load('C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/CP11_T3_Opt7_wkR_0.05.mat');
v9 = load('C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/CP11_T3_Opt7_wkR_0.005.mat');
v10 = load('C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/CP11_T3_Opt7_wkR_0.7.mat');
v11 = load('C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/CP11_T3_Opt7_wkR_0.07.mat');
v12 = load('C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/CP11_T3_Opt7_wkR_0.007.mat');
v13 = load('C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/CP11_T3_Opt7_wkR_1.mat');
v14 = load('C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/CP11_T3_Opt7_wkR_10.mat');

gain     = [v1.R.kR v2.R.kR v3.R.kR v4.R.kR v5.R.kR v6.R.kR v7.R.kR v8.R.kR v9.R.kR v10.R.kR v11.R.kR v12.R.kR v13.R.kR v14.R.kR ];
kost     = [v1.R.J v2.R.J v3.R.J v4.R.J v5.R.J v6.R.J v7.R.J v8.R.J v9.R.J v10.R.J v11.R.J v12.R.J v13.R.J v14.R.J ];
kost_no  = kost - gain.*[0.1 0.01 0.001 0.2 0.02 0.002 0.5 0.05 0.005 0.7 0.07 0.007 1 10]; 

figure()
scatter(kost_no,gain,'Filled'); hold on
box off; title('Pareto CP11 T3'); ylabel('Gain'); xlabel('Kost without gain'); 
saveas(gcf,'C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/Figures/CP11_T3_Pareto.png');
saveas(gcf,'C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/Figures/CP11_T3_Pareto.fig');

figure()
scatter(kost_no,(gain).^2,'Filled'); hold on
box off; title('Pareto CP11 T3'); ylabel('Gain'); xlabel('Kost without gain'); 
saveas(gcf,'C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/Figures/CP11_T3_Pareto2.png');
saveas(gcf,'C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/Figures/CP11_T3_Pareto2.fig');
