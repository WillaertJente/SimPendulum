%% Script to compare simulated and experimental reflex activity of the quadriceps muscle 
% Jente Willaert - 22 april 2021
clear all; close all; clc
%% Load data
map_traj = 'C:\Users\u0125183\Box\PhD 1\Simulations Pendulum Test\Results/'; 
map_emg  = 'C:\Users\u0125183\Box\PhD 1\Pendulum test\Data\CP16\Pend_20190729_/';
map_out  = 'C:\Users\u0125183\Box\PhD 1\Simulations Pendulum Test\Spiermodel\Plots/';
subject  = 'CP16'; 

for i = 1:20
    FullFileName = [map_traj,'Result_',subject,'_T',num2str(i),'_Fmv.mat']
    if exist(FullFileName)
        % Load simulated and experimental trajectories
        data_traj = load(FullFileName); 
        
        % Load experimental EMG data 
        EMGFileName = [map_emg,'/EMG/SynctEMG_v2/SynctEMG_trial',num2str(i),'.xls']; 
        data_emg    = xlsread(EMGFileName,'1000Hz');
        
        start_emg   = find(data_emg(:,1)==data_traj.tvect(1))
        
        % Calculate simulated reflex activity
        Fsrs_d_calc = zeros(length(data_traj.q_exp),1); 
        for s = 1:length(data_traj.q_exp)
            Fsrs_d_calc(s) = (data_traj.sol_a_ext(s)-data_traj.sol_aext0)/data_traj.sol_Rk;
        end
        Reflex_calc = Fsrs_d_calc*data_traj.sol_Rk;
        
%% Plot
        figure(i)
        suptitle([subject, ' T', num2str(i)])
        subplot(311)
        plot(data_traj.tvect,data_traj.q_exp*180/pi,'k','LineWidth',1.5); hold on;
        plot(data_traj.tvect,data_traj.sol_x*180/pi,'r','LineWidth',1.5); hold on; 
        legend({'Exp','Sim'}); ylabel('[{\circ}]'); box off
        
        subplot(312)
        plot(data_emg(start_emg:end,1),data_emg(start_emg:end,2),'LineWidth',1.5); hold on; 
        plot(data_emg(start_emg:end,1),data_emg(start_emg:end,3),'LineWidth',1.5); hold on;
        legend({'RF','VL'}); ylabel('[V]'); box off; title('Experimental EMG')
        
        subplot(313)
        plot(data_traj.tvect,Reflex_calc,'k','LineWidth',1.5); hold on; 
        title('Simulated reflex activity'); ylabel([]); box off 
        legend({['Rk=',num2str(data_traj.sol_Rk)]})
    
        saveas(gcf,[map_out,subject,'/EMG_T',num2str(i),'.png'])
    end
end