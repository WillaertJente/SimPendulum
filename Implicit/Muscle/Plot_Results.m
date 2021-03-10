%% Script to plot results of simulations
% Jente Willaert
% 9/03/2021
clear all
close all
clc

% Import data
subj = 'CP15';
map = 'C:\Users\u0125183\Box\PhD 1\Simulations Pendulum Test\Results/';

for i = 1:20
    Name = [map,'Result_',char(subj),'_T',num2str(i),'_Refl4.mat'];
    if exist(Name)
        load(Name)
        
        Total = sol_FT_ext.*sol_ma_ext + sol_FT_flex.*sol_ma_flex;
        tvect = 0.005:0.005:length(q_exp)/200;
        figure()
        suptitle(['Trial ',num2str(i)])
        
        subplot(611)
        plot(tvect,q_exp*180/pi,'k','LineWidth',1.5); hold on
        plot(tvect,sol_x*180/pi,'LineWidth',1.5); hold on
        ylabel('KA ({\circ})'); box off; legend('Exp','Sim')
        
        subplot(612)
        plot(tvect,sol_FT_ext,'LineWidth',1.5); hold on;
        plot(tvect,sol_FT_flex,'LineWidth',1.5); hold on;
        ylabel('FT (N)'); box off; legend('RF','BF');
        
        subplot(613)
        plot(tvect,sol_Fpe_ext,'LineWidth',1.5); hold on;
        plot(tvect,sol_Fpe_flex,'LineWidth',1.5); hold on;
        ylabel('Fpe (Norm)'); box off; legend('RF','BF');
        
        subplot(614)
        plot(tvect,sol_Fsrs,'LineWidth',1.5); hold on
        ylabel('Fsrs (Norm)'); box off;
        
        subplot(615)
        plot(tvect,sol_a_ext,'LineWidth',1.5); hold on;
        plot(tvect,sol_a_flex*ones(1,length(tvect)),'LineWidth',1.5); hold on
        ylabel('a'); box off; legend('RF','BF');
        
        subplot(616)
        plot(tvect,sol_act,'Color',[1 0 0 ],'LineWidth',1.5); hold on;
        plot(tvect,Total,'Color','k','LineWidth',1.5); hold on;
        ylabel('Nm'); box off; legend('Actuator','Muscle torque');
        
        saveas(gcf,['C:\Users\u0125183\Box\PhD 1\Simulations Pendulum Test/Spiermodel/Plots/',char(subj),'/char(subj)_T',num2str(i),'_Refl4.png'])
    end
    
end

%% Reflex vs. no Reflex
subj = 'CP15';
map = 'C:\Users\u0125183\Box\PhD 1\Simulations Pendulum Test\Results/';

for i = [4]
    Name_r = [map,'Result_',char(subj),'_T',num2str(i),'_Refl4.mat'];
    Name_nr= [map,'Result_',char(subj),'_T',num2str(i),'_Refl1.mat'];
    if exist(Name_nr)
        
        load(Name_nr)
        tvect = 0.005:0.005:length(q_exp)/200;
        figure()
        plot(tvect,q_exp*180/pi, 'k', 'LineWidth', 1.5); hold on;
        plot(tvect,sol_x*180/pi, 'LineWidth', 1.5); hold on
        
        load(Name_r)
        plot(tvect,sol_x*180/pi,'LineWidth',1.5); hold on;
        title(['Trial ',num2str(i)]); box off; ylabel('KA ({\circ})');
        legend({'Exp','No Reflex', 'Reflex'})
        
        saveas(gcf,['C:\Users\u0125183\Box\PhD 1\Simulations Pendulum Test/Spiermodel/Plots/',char(subj),'/',char(subj),'_T',num2str(i),'_RvsNR.png'])
    else
        disp(['Trial', num2str(i),'niet geconvergeerd'])
    end
    
end
