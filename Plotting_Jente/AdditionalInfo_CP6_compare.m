%% Script to plot additional trajectories 

% Trial 1
a = load('C:\Users\u0125183\Box\PhD 1\Simulations Pendulum Test\Results/Result_CP6_T2_FMv.mat'); 
b = load('C:\Users\u0125183\Box\PhD 1\Simulations Pendulum Test\Results/Result_CP6_T2_ScaledTorso2.mat'); 
t = a.tvect-a.tvect(1); 

figure(1)
subplot(221)
plot(t,a.q_exp*180/pi,'k','LineWidth',1.5); hold on; 
plot(t,a.sol_x*180/pi,'Color',[2 118 161]./255,'LineWidth',1.5); hold on;
plot(t,a.sol_x*180/pi,'Color',[52 186 235]./255,'LineWidth',1.5); hold on;
ylabel('Knee angle ({\circ})'); xlabel('Time (s)'); legend({'Exp','Fmv','ScaledTorso'})
title('Trial 1'); box off
axis([0 13 -150 0])

figure(2)
subplot(231)
bar(1,a.sol_aext0,0.5,'EdgeColor',[2 118 161]./255,'FaceColor',[2 118 161]./255); hold on; 
bar(2,b.sol_aext0,0.5,'EdgeColor',[52 186 235]./255,'FaceColor',[52 186 235]./255); hold on
bar(4,a.sol_a_flex, 0.5,'EdgeColor',[2 118 161]./255,'FaceColor',[2 118 161]./255); hold on; 
bar(5,b.sol_a_flex,0.5,'EdgeColor',[52 186 235]./255,'FaceColor',[52 186 235]./255); hold on
xticks([1.5 4.5])
xticklabels({'Ext','Flexor'})
title('Baseline muscle tone'); box off; 
axis([0 6 0 0.02])

subplot(232)
bar(1,0.15-a.sol_kFpe_ext,0.5,'EdgeColor',[2 118 161]./255,'FaceColor',[2 118 161]./255); hold on; 
bar(2,0.15-b.sol_kFpe_ext,0.5,'EdgeColor',[52 186 235]./255,'FaceColor',[52 186 235]./255); hold on; 
bar(4,0.15-a.sol_kFpe_flex, 0.5,'EdgeColor',[2 118 161]./255,'FaceColor',[2 118 161]./255); hold on; 
bar(5,0.15-b.sol_kFpe_flex, 0.5,'EdgeColor',[52 186 235]./255,'FaceColor',[52 186 235]./255); hold on; 
xticks([1.5 4.5])
xticklabels({'Ext','Flexor'})
title('Passive muscle stiffness');  box off;
axis([0 6 0 0.1])

subplot(233)
bar(1,a.sol_Rk,0.5,'EdgeColor',[2 118 161]./255,'FaceColor',[2 118 161]./255); hold on;
bar(2,b.sol_Rk,0.5,'EdgeColor',[52 186 235]./255,'FaceColor',[2 118 161]./255); hold on;
title('Reflex gain'); box off
xticks([1])
xticklabels({'Ext'})
axis([0 3 0 4])

% Trial 2
a = load('C:\Users\u0125183\Box\PhD 1\Simulations Pendulum Test\Results/Result_CP6_T3_FMv.mat');
b = load('C:\Users\u0125183\Box\PhD 1\Simulations Pendulum Test\Results/Result_CP6_T3_ScaledTorso2.mat');
t = a.tvect-a.tvect(1); 

figure(3)
subplot(221)
plot(t,a.q_exp*180/pi,'k','LineWidth',1.5); hold on; 
plot(t,a.sol_x*180/pi,'Color',[2 118 161]./255,'LineWidth',1.5); hold on;
plot(t,b.sol_x*180/pi,'Color',[52 186 235]./255,'LineWidth',1.5); hold on;
ylabel('Knee angle ({\circ})'); xlabel('Time (s)'); legend({'Exp','Fmv','ScaledTorso'})
title('Trial 2'); box off
axis([0 13 -150 0])

figure(4)
subplot(231)
bar(1,a.sol_aext0,0.5,'EdgeColor',[2 118 161]./255,'FaceColor',[2 118 161]./255); hold on; 
bar(2,b.sol_aext0,0.5,'EdgeColor',[52 186 235]./255,'FaceColor',[52 186 235]./255); hold on; 
bar(4,a.sol_a_flex, 0.5,'EdgeColor',[2 118 161]./255,'FaceColor',[2 118 161]./255); hold on; 
bar(5,b.sol_a_flex, 0.5,'EdgeColor',[52 186 235]./255,'FaceColor',[52 186 235]./255); hold on; 
xticks([1.5 4.5])
xticklabels({'Ext','Flexor'})
title('Baseline muscle tone');  box off; 
axis([0 6 0 0.02])

subplot(232)
bar(1,0.15-a.sol_kFpe_ext,0.5,'EdgeColor',[2 118 161]./255,'FaceColor',[2 118 161]./255); hold on; 
bar(2,0.15-b.sol_kFpe_ext,0.5,'EdgeColor',[52 186 235]./255,'FaceColor',[52 186 235]./255); hold on; 
bar(4,0.15-a.sol_kFpe_flex, 0.5,'EdgeColor',[2 118 161]./255,'FaceColor',[2 118 161]./255); hold on; 
bar(5,0.15-b.sol_kFpe_flex, 0.5,'EdgeColor',[52 186 235]./255,'FaceColor',[52 186 235]./255); hold on; 
xticks([1.5 4.5])
xticklabels({'Ext','Flexor'})
title('Passive muscle stiffness');  box off;
axis([0 6 0 0.1])

subplot(233)
bar(1,a.sol_Rk,0.5,'EdgeColor',[2 118 161]./255,'FaceColor',[2 118 161]./255); hold on; 
bar(2,a.sol_Rk,0.5,'EdgeColor',[52 186 235]./255,'FaceColor',[52 186 235]./255); hold on; 
title('Reflex gain'); box off
xticks([1.5])
xticklabels({'Ext'})
axis([0 3 0 4])

% Trial 3
a = load('C:\Users\u0125183\Box\PhD 1\Simulations Pendulum Test\Results/Result_CP6_T4_FMv.mat'); 
b = load('C:\Users\u0125183\Box\PhD 1\Simulations Pendulum Test\Results/Result_CP6_T4_ScaledTorso.mat'); 
t = a.tvect-a.tvect(1); 

figure(5)
subplot(221)
plot(t,a.q_exp*180/pi,'k','LineWidth',1.5); hold on; 
plot(t,a.sol_x*180/pi,'Color',[2 118 161]./255,'LineWidth',1.5); hold on;
plot(t,b.sol_x*180/pi,'Color',[52 186 235]./255,'LineWidth',1.5); hold on;
ylabel('Knee angle ({\circ})'); xlabel('Time (s)'); legend({'Exp','Sim'})
title('Trial 3'); box off
axis([0 13 -150 0])

figure(6)
subplot(231)
bar(1,a.sol_aext0,0.5,'EdgeColor',[2 118 161]./255,'FaceColor',[2 118 161]./255); hold on; 
bar(2,b.sol_aext0,0.5,'EdgeColor',[52 186 235]./255,'FaceColor',[52 186 235]./255); hold on; 
bar(4,a.sol_a_flex, 0.5,'EdgeColor',[2 118 161]./255,'FaceColor',[2 118 161]./255); hold on; 
bar(5,b.sol_a_flex, 0.5,'EdgeColor',[52 186 235]./255,'FaceColor',[52 186 235]./255); hold on; 
xticks([1.5 4.5])
xticklabels({'Ext','Flexor'})
title('Baseline muscle tone');  box off; 
axis([0 6 0 0.02])

subplot(232)
bar(1,0.15-a.sol_kFpe_ext,0.5,'EdgeColor',[2 118 161]./255,'FaceColor',[2 118 161]./255); hold on; 
bar(2,0.15-b.sol_kFpe_ext,0.5,'EdgeColor',[52 186 235]./255,'FaceColor',[52 186 235]./255); hold on; 
bar(4,0.15-a.sol_kFpe_flex, 0.5,'EdgeColor',[2 118 161]./255,'FaceColor',[2 118 161]./255); hold on; 
bar(5,0.15-b.sol_kFpe_flex, 0.5,'EdgeColor',[52 186 235]./255,'FaceColor',[52 186 235]./255); hold on; 
xticks([1.5 4.5])
xticklabels({'Ext','Flexor'})
title('Passive muscle stiffness');  box off;
axis([0 6 0 0.1])

subplot(233)
bar(1,a.sol_Rk,0.5,'EdgeColor',[2 118 161]./255,'FaceColor',[2 118 161]./255); hold on; 
bar(2,a.sol_Rk,0.5,'EdgeColor',[52 186 235]./255,'FaceColor',[52 186 235]./255); hold on; 
title('Reflex gain'); box off
xticks([1.5])
xticklabels({'Ext'})
axis([0 3 0 4])

% Trial 2
a = load('C:\Users\u0125183\Box\PhD 1\Simulations Pendulum Test\Results/Result_CP6_T6_FMv.mat')
b = load('C:\Users\u0125183\Box\PhD 1\Simulations Pendulum Test\Results/Result_CP6_T6_ScaledTorso.mat')
t = a.tvect-a.tvect(1); 

figure(7)
subplot(221)
plot(t,a.q_exp*180/pi,'k','LineWidth',1.5); hold on; 
plot(t,a.sol_x*180/pi,'Color',[2 118 161]./255,'LineWidth',1.5); hold on;
plot(t,b.sol_x*180/pi,'Color',[52 186 235]./255,'LineWidth',1.5); hold on;
ylabel('Knee angle ({\circ})'); xlabel('Time (s)'); legend({'Exp','Sim'})
title('Trial 4'); box off
axis([0 13 -150 0])

figure(8)
subplot(231)
bar(1,a.sol_aext0,0.5,'EdgeColor',[2 118 161]./255,'FaceColor',[2 118 161]./255); hold on; 
bar(2,b.sol_aext0,0.5,'EdgeColor',[52 186 235]./255,'FaceColor',[52 186 235]./255); hold on; 
bar(4,a.sol_a_flex, 0.5,'EdgeColor',[2 118 161]./255,'FaceColor',[2 118 161]./255); hold on; 
bar(5,b.sol_a_flex, 0.5,'EdgeColor',[52 186 235]./255,'FaceColor',[52 186 235]./255); hold on; 
xticks([1.5 4.5])
xticklabels({'Ext','Flexor'})
title('Baseline muscle tone');  box off; 
axis([0 6 0 0.02])

subplot(232)
bar(1,0.15-a.sol_kFpe_ext,0.5,'EdgeColor',[2 118 161]./255,'FaceColor',[2 118 161]./255); hold on; 
bar(2,0.15-b.sol_kFpe_ext,0.5,'EdgeColor',[52 186 235]./255,'FaceColor',[52 186 235]./255); hold on; 
bar(4,0.15-a.sol_kFpe_flex, 0.5,'EdgeColor',[2 118 161]./255,'FaceColor',[2 118 161]./255); hold on; 
bar(5,0.15-b.sol_kFpe_flex, 0.5,'EdgeColor',[52 186 235]./255,'FaceColor',[52 186 235]./255); hold on; 
xticks([1.5 4.5])
xticklabels({'Ext','Flexor'})
title('Passive muscle stiffness');  box off;
axis([0 6 0 0.1])

subplot(233)
bar(1,a.sol_Rk,0.5,'EdgeColor',[2 118 161]./255,'FaceColor',[2 118 161]./255); hold on; 
bar(2,a.sol_Rk,0.5,'EdgeColor',[52 186 235]./255,'FaceColor',[52 186 235]./255); hold on; 
title('Reflex gain'); box off
xticks([1.5])
xticklabels({'Ext'})
axis([0 3 0 4])