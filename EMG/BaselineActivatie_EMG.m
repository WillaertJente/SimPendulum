%% BaselineActivatie EMG 
clear all; close all; clc

% Input
for t = [1 2]
    subj  = 'TD5';
    file_in   = ['C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\EMG\',char(subj),'/T',num2str(t),'.xlsx'];
    file_out  = ['C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\EMG\',char(subj),'/T',num2str(t),'_Sim.xlsx'];

    all   = xlsread(file_in);
    sim   = xlsread(file_out);
    start = sim(1,1)*1000;

    plot(all(start-500:start+2000)); hold on
    line([500 500],[0 0.00001],'Color',[0.7 0.7 0.7],'LineWidth',1.5)
    
    avg(t) = mean(all(start-500:start)); 
    clear start
end

%xlswrite(['C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\EMG/',char(subj),'/BaselineActivation_prestretch.xlsx'],avg)



