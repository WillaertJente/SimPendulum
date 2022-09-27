% TD5 T1 
td5 = load('C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/TD5_T1_Opt12_IG_set5_Bounds_Set1.mat');

figure()
subplot(221)
plot(td5.R.exp.qspline,'k','LineWidth',1.5); hold on
plot(td5.R.x,'r','LineWidth',1.2); hold on
title('Healthy example'); box off
subplot(256,'replace')
bar(1,td5.R.a_ext); hold on
box off; title('a extensor'); ylim([0 0.02])
subplot(257)
bar(1,td5.R.a_flex); hold on
box off; title('a flexor'); ylim([0 0.02])
subplot(258)
bar(1,td5.R.kR); hold on
box off; title('Reflex gain'); ylim([0 5])
subplot(259)
bar(1,td5.R.kFpe); hold on
box off; title('Passive stiffness'); ylim([0 0.22])
subplot(2,5,10)
bar(1,td5.R.B); hold on
box off; title('Damping'); ylim([0 0.1])

% CP9 T5 
cp9 = load('C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/CP9_T5_Opt12_IG_set4_Bounds_Set1.mat');

figure()
subplot(221)
plot(cp9.R.exp.qspline,'k','LineWidth',1.5); hold on
plot(cp9.R.x,'r','LineWidth',1.2); hold on
title('No spasticity'); box off
subplot(256,'replace')
bar(1,cp9.R.a_ext); hold on
box off; title('a extensor'); ylim([0 0.02])
subplot(257)
bar(1,cp9.R.a_flex); hold on
box off; title('a flexor'); ylim([0 0.02])
subplot(258)
bar(1,cp9.R.kR); hold on
box off; title('Reflex gain'); ylim([0 5])
subplot(259)
bar(1,cp9.R.kFpe); hold on
box off; title('Passive stiffness'); ylim([0 0.22])
subplot(2,5,10)
bar(1,cp9.R.B); hold on
box off; title('Damping'); ylim([0 0.1])


% CP8 T2 
cp8 = load('C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/CP8_T2_Opt12_IG_set7_Bounds_Set1.mat');


figure()
subplot(221)
plot(cp8.R.exp.qspline,'k','LineWidth',1.5); hold on
plot(cp8.R.x,'r','LineWidth',1.2); hold on
title('High reflex actiivy'); box off
subplot(256,'replace')
bar(1,cp8.R.a_ext); hold on
box off; title('a extensor'); ylim([0 0.02])
subplot(257)
bar(1,cp8.R.a_flex); hold on
box off; title('a flexor'); ylim([0 0.02])
subplot(258)
bar(1,cp8.R.kR); hold on
box off; title('Reflex gain'); ylim([0 5])
subplot(259)
bar(1,cp8.R.kFpe); hold on
box off; title('Passive stiffness'); ylim([0 0.22])
subplot(2,5,10)
bar(1,cp8.R.B); hold on
box off; title('Damping'); ylim([0 0.1])


% CP11 T3 
cp11 = load('C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/CP11_T3_Opt12_IG_set8_Bounds_Set1.mat');


figure()
subplot(221)
plot(cp11.R.exp.qspline,'k','LineWidth',1.5); hold on
plot(cp11.R.x,'r','LineWidth',1.2); hold on
title('Neural and non-neural'); box off
subplot(256,'replace')
bar(1,cp11.R.a_ext); hold on
box off; title('a extensor'); ylim([0 0.02])
subplot(257)
bar(1,cp11.R.a_flex); hold on
box off; title('a flexor'); ylim([0 0.02])
subplot(258)
bar(1,cp11.R.kR); hold on
box off; title('Reflex gain'); ylim([0 5])
subplot(259)
bar(1,cp11.R.kFpe); hold on
box off; title('Passive stiffness'); ylim([0 0.22])
subplot(2,5,10)
bar(1,cp11.R.B); hold on
box off; title('Damping'); ylim([0 0.1])