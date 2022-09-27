%% BaselineActivatie EMG 
clear all; close all; clc

% Input
for t = [8 9 10 11 ]
    subj  = 'CP11';
    file_in   = ['C:\Users\u0125183\OneDrive - KU Leuven\Pendulum test\Data\',char(subj),'/EMG/SynctEMG_trial',num2str(t),'.xls'];
    
    all   = xlsread(file_in,'Original');
    
    plot(all(:,1),all(:,2)); hold on
    box off; title('baseline EMG - PM')
    ylim([-0.00001 0.00007])
    clear start
end

%xlswrite(['C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\EMG/',char(subj),'/BaselineActivation_prestretch.xlsx'],avg)



