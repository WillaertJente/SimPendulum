%% Plot results
% CP 2
t4 = load('C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/CP2_T4_Opt12_IG_set4_Bounds_Set2.mat');
t5 = load('C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/CP2_T5_Opt12_IG_set4_Bounds_Set3.mat');
t6 = load('C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/CP2_T6_Opt12_IG_set4_Bounds_Set3.mat');
t7 = load('C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/CP2_T7_Opt12_IG_set4_Bounds_Set2.mat');

figure()
subplot(251)
plot(t4.R.exp.qspline,'k'); hold on
plot(t4.R.x); hold on
box off; title('CP2 - T4')

subplot(252)
plot(t5.R.exp.qspline,'k'); hold on
plot(t5.R.x); hold on
box off; title('CP2 - T5')

subplot(253)
plot(t6.R.exp.qspline,'k'); hold on
plot(t6.R.x); hold on
box off; title('CP2 - T6')

subplot(254)
plot(t7.R.exp.qspline,'k'); hold on
plot(t7.R.x); hold on
box off; title('CP2 - T7')

subplot(256)
bar(1,t4.R.a_ext); hold on
bar(2,t5.R.a_ext); hold on
bar(3,t6.R.a_ext); hold on
bar(4,t7.R.a_ext); hold on
box off; title('a ext')

subplot(257)
bar(1,t4.R.a_flex); hold on
bar(2,t5.R.a_flex); hold on
bar(3,t6.R.a_flex); hold on
bar(4,t7.R.a_flex); hold on
box off; title('a flex')

subplot(258)
bar(1,t4.R.kR); hold on
bar(2,t5.R.kR); hold on
bar(3,t6.R.kR); hold on
bar(4,t7.R.kR); hold on
box off; title('R')

subplot(259)
bar(1,t4.R.kFpe); hold on
bar(2,t5.R.kFpe); hold on
bar(3,t6.R.kFpe); hold on
bar(4,t7.R.kFpe); hold on
box off; title('kFpe')

subplot(2,5,10)
bar(1,t4.R.B); hold on
bar(2,t5.R.B); hold on
bar(3,t6.R.B); hold on
bar(4,t7.R.B); hold on
box off; title('B')

% CP 
t3 = load('C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/CP4_T3_Opt12_IG_set8_Bounds_Set3.mat');
t4 = load('C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/CP4_T4_Opt12_IG_set8_Bounds_Set3.mat');

figure()
subplot(251)
plot(t3.R.exp.qspline,'k'); hold on
plot(t3.R.x); hold on
box off; title('CP4 - T3')

subplot(252)
plot(t4.R.exp.qspline,'k'); hold on
plot(t4.R.x); hold on
box off; title('CP4 - T4')

subplot(256)
bar(1,t3.R.a_ext); hold on
bar(2,t4.R.a_ext); hold on
box off; title('a ext')

subplot(257)
bar(1,t3.R.a_flex); hold on
bar(2,t4.R.a_flex); hold on
box off; title('a flex')

subplot(258)
bar(1,t3.R.kR); hold on
bar(2,t4.R.kR); hold on
box off; title('R')

subplot(259)
bar(1,t3.R.kFpe); hold on
bar(2,t4.R.kFpe); hold on
box off; title('kFpe')

subplot(2,5,10)
bar(1,t3.R.B); hold on
bar(2,t4.R.B); hold on
box off; title('B')

% CP 8
t2 = load('C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/CP8_T2_Opt12_IG_set7_Bounds_Set1.mat');
t3 = load('C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/CP8_T3_Opt12_IG_set7_Bounds_Set1.mat');
t4 = load('C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/CP8_T4_Opt12_IG_set7_Bounds_Set1.mat');
t5 = load('C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/CP8_T5_Opt12_IG_set7_Bounds_Set1.mat');

figure()
subplot(251)
plot(t2.R.exp.qspline,'k'); hold on
plot(t2.R.x); hold on
box off; title('CP8 - T2')

subplot(252)
plot(t3.R.exp.qspline,'k'); hold on
plot(t3.R.x); hold on
box off; title('CP8 - T3')

subplot(253)
plot(t4.R.exp.qspline,'k'); hold on
plot(t4.R.x); hold on
box off; title('CP8 - T4')

subplot(254)
plot(t5.R.exp.qspline,'k'); hold on
plot(t5.R.x); hold on
box off; title('CP8 - T5')

subplot(256)
bar(1,t2.R.a_ext); hold on
bar(2,t3.R.a_ext); hold on
bar(3,t4.R.a_ext); hold on
bar(4,t5.R.a_ext); hold on
box off; title('a ext')

subplot(257)
bar(1,t2.R.a_flex); hold on
bar(2,t3.R.a_flex); hold on
bar(3,t4.R.a_flex); hold on
bar(4,t5.R.a_flex); hold on
box off; title('a flex')

subplot(258)
bar(1,t2.R.kR); hold on
bar(2,t3.R.kR); hold on
bar(3,t4.R.kR); hold on
bar(4,t5.R.kR); hold on
box off; title('R')

subplot(259)
bar(1,t2.R.kFpe); hold on
bar(2,t3.R.kFpe); hold on
bar(3,t4.R.kFpe); hold on
bar(4,t5.R.kFpe); hold on
box off; title('kFpe')

subplot(2,5,10)
bar(1,t2.R.B); hold on
bar(2,t3.R.B); hold on
bar(3,t4.R.B); hold on
bar(4,t5.R.B); hold on
box off; title('B')

% CP 9
t2 = load('C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/CP9_T2_Opt12_IG_set8_Bounds_Set2.mat');
t3 = load('C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/CP9_T3_Opt12_IG_set8_Bounds_Set2.mat');
t4 = load('C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/CP9_T4_Opt12_IG_set8_Bounds_Set3.mat');
t5 = load('C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/CP9_T5_Opt12_IG_set4_Bounds_Set1.mat');
t6 = load('C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/CP9_T6_Opt12_IG_set8_Bounds_Set2.mat');

figure()
subplot(251)
plot(t2.R.exp.qspline,'k'); hold on
plot(t2.R.x); hold on
box off; title('CP9 - T2')

subplot(252)
plot(t3.R.exp.qspline,'k'); hold on
plot(t3.R.x); hold on
box off; title('CP9 - T3')

subplot(253)
plot(t4.R.exp.qspline,'k'); hold on
plot(t4.R.x); hold on
box off; title('CP9 - T4')

subplot(254)
plot(t5.R.exp.qspline,'k'); hold on
plot(t5.R.x); hold on
box off; title('CP9 - T5')

subplot(255)
plot(t6.R.exp.qspline,'k'); hold on
plot(t6.R.x); hold on
box off; title('CP9 - T6')

subplot(256)
bar(1,t2.R.a_ext); hold on
bar(2,t3.R.a_ext); hold on
bar(3,t4.R.a_ext); hold on
bar(4,t5.R.a_ext); hold on
bar(5,t6.R.a_ext); hold on
box off; title('a ext')

subplot(257)
bar(1,t2.R.a_flex); hold on
bar(2,t3.R.a_flex); hold on
bar(3,t4.R.a_flex); hold on
bar(4,t5.R.a_flex); hold on
bar(5,t6.R.a_flex); hold on
box off; title('a flex')

subplot(258)
bar(1,t2.R.kR); hold on
bar(2,t3.R.kR); hold on
bar(3,t4.R.kR); hold on
bar(4,t5.R.kR); hold on
bar(5,t6.R.kR); hold on
box off; title('R')

subplot(259)
bar(1,t2.R.kFpe); hold on
bar(2,t3.R.kFpe); hold on
bar(3,t4.R.kFpe); hold on
bar(4,t5.R.kFpe); hold on
bar(5,t6.R.kFpe); hold on
box off; title('kFpe')

subplot(2,5,10)
bar(1,t2.R.B); hold on
bar(2,t3.R.B); hold on
bar(3,t4.R.B); hold on
bar(4,t5.R.B); hold on
bar(4,t6.R.B); hold on
box off; title('B')

% CP 10
t12 = load('C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/CP10_T12_Opt12_IG_set6_Bounds_Set2.mat');
t13 = load('C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/CP10_T13_Opt12_IG_set6_Bounds_Set2.mat');
t14 = load('C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/CP10_T14_Opt12_IG_set9_Bounds_Set3.mat');

figure()
subplot(251)
plot(t12.R.exp.qspline,'k'); hold on
plot(t12.R.x); hold on
box off; title('CP10 - T12')

subplot(252)
plot(t13.R.exp.qspline,'k'); hold on
plot(t13.R.x); hold on
box off; title('CP10 - T13')

subplot(253)
plot(t14.R.exp.qspline,'k'); hold on
plot(t14.R.x); hold on
box off; title('CP10 - T14')

subplot(256)
bar(1,t12.R.a_ext); hold on
bar(2,t13.R.a_ext); hold on
bar(3,t14.R.a_ext); hold on
box off; title('a ext')

subplot(257)
bar(1,t12.R.a_flex); hold on
bar(2,t13.R.a_flex); hold on
bar(3,t14.R.a_flex); hold on
box off; title('a flex')

subplot(258)
bar(1,t12.R.kR); hold on
bar(2,t13.R.kR); hold on
bar(3,t14.R.kR); hold on
box off; title('R')

subplot(259)
bar(1,t12.R.kFpe); hold on
bar(2,t13.R.kFpe); hold on
bar(3,t14.R.kFpe); hold on
box off; title('kFpe')

subplot(2,5,10)
bar(1,t12.R.B); hold on
bar(2,t13.R.B); hold on
bar(3,t14.R.B); hold on
box off; title('B')

% CP 11
t2 = load('C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/CP11_T2_Opt12_IG_set8_Bounds_Set3.mat');
t3 = load('C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/CP11_T3_Opt12_IG_set8_Bounds_Set1.mat');
t4 = load('C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/CP11_T4_Opt12_IG_set8_Bounds_Set1.mat');

figure()
subplot(251)
plot(t2.R.exp.qspline,'k'); hold on
plot(t2.R.x); hold on
box off; title('CP11 - T2')

subplot(252)
plot(t3.R.exp.qspline,'k'); hold on
plot(t3.R.x); hold on
box off; title('CP11 - T3')

subplot(253)
plot(t4.R.exp.qspline,'k'); hold on
plot(t4.R.x); hold on
box off; title('CP11 - T4')

subplot(256)
bar(1,t2.R.a_ext); hold on
bar(2,t3.R.a_ext); hold on
bar(3,t4.R.a_ext); hold on
box off; title('a ext')

subplot(257)
bar(1,t2.R.a_flex); hold on
bar(2,t3.R.a_flex); hold on
bar(3,t4.R.a_flex); hold on
box off; title('a flex')

subplot(258)
bar(1,t2.R.kR); hold on
bar(2,t3.R.kR); hold on
bar(3,t4.R.kR); hold on
box off; title('R')

subplot(259)
bar(1,t2.R.kFpe); hold on
bar(2,t3.R.kFpe); hold on
bar(3,t4.R.kFpe); hold on
box off; title('kFpe')

subplot(2,5,10)
bar(1,t2.R.B); hold on
bar(2,t3.R.B); hold on
bar(3,t4.R.B); hold on
box off; title('B')

% TD 5
t1 = load('C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/TD5_T1_Opt12_IG_set5_Bounds_Set1.mat');
t2 = load('C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/TD5_T2_Opt12_IG_set5_Bounds_Set1.mat');

figure()
subplot(251)
plot(t1.R.exp.qspline,'k'); hold on
plot(t1.R.x); hold on
box off; title('TD5 - T1')

subplot(252)
plot(t2.R.exp.qspline,'k'); hold on
plot(t2.R.x); hold on
box off; title('TD5 - T2')

subplot(256)
bar(1,t1.R.a_ext); hold on
bar(2,t2.R.a_ext); hold on
box off; title('a ext')

subplot(257)
bar(1,t1.R.a_flex); hold on
bar(2,t2.R.a_flex); hold on
box off; title('a flex')

subplot(258)
bar(1,t1.R.kR); hold on
bar(2,t2.R.kR); hold on
box off; title('R')

subplot(259)
bar(1,t1.R.kFpe); hold on
bar(2,t2.R.kFpe); hold on
box off; title('kFpe')

subplot(2,5,10)
bar(1,t1.R.B); hold on
bar(2,t2.R.B); hold on
box off; title('B')