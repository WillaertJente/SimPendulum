%% BaselineActivatie EMG 
clear all; close all; clc

% Input
for t = [3 4]
    subj  = 'CP4';
    file_in   = ['C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\EMG\',char(subj),'/T',num2str(t),'.xlsx'];
    file_out  = ['C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\EMG\',char(subj),'/T',num2str(t),'_Sim.xlsx'];

    all   = xlsread(file_in);
    sim   = xlsread(file_out);
    start = sim(1,1)*1000;

    plot(all(:,1)); hold on
    %line([start start],[0 0.000005],'LineWidth',2); hold on
    box off; title('baseline EMG')
    ylim([-0.000001 0.000006])
    clear start
end

%xlswrite(['C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\EMG/',char(subj),'/BaselineActivation_prestretch.xlsx'],avg)



