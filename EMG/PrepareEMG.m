%% Prepare EMG 
clear all; close all; clc
% Input
t     = 4; 
subj  = 'CP11';
file_in   = ['C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\EMG\',char(subj),'/T',num2str(t),'.xlsx'];
file_out  = ['C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\EMG\',char(subj),'/T',num2str(t),'_Sim.xlsx'];
file_bk   = ['C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Implicit\Muscle\Experimental data\', char(subj),'/BK_Trial',num2str(t),'.mat'];


% Read data
data_emg   = xlsread(file_in); 
data_bk    = load(file_bk); 
tspan      = [data_bk.data(1,1) data_bk.data(end,1)];

% Selected data
sel        = data_emg(tspan(1)*1000:tspan(2)*1000,1:5);
time       = tspan(1):0.001:tspan(2);

% Write to excel 
Headers = {'Time','RF','VL','VM','BF','ST'}; 

xlswrite(file_out,Headers,'Sheet1','A1');
xlswrite(file_out,time','Sheet1','A2');
xlswrite(file_out,sel,'Sheet1','B2')
