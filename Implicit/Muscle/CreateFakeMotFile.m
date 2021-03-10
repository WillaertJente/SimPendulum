%% Create fake mot file for OpenSim

%% Hip at mean of first 100 frames
addpath('C:\Users\u0125183\Documents\MATLAB\Sam/WriteMotFile')
leg = 0; % 0 is rechts, 1 is links
% Read existing mot file
[header, names, data, fpath] = SIMM_ReadMotion('C:\Users\u0125183\Box\PhD 1\Simulations Pendulum Test\Data\TD12/IK_Trial2.mot');

% Read hip flexion
col_hipf_r = find(strcmp('hip_flexion_r',names))
col_hipf_l = find(strcmp('hip_flexion_l',names))

if leg == 0;
    figure()
    plot(data(:,col_hipf_r))
    mean_hip_flexion = mean(data(1:100,col_hipf_r))
else
    figure()
    plot(data(:,col_hipf_l))
    mean_hip_flexion = mean(data(1:100,col_hipf_l))
end

% new data 
knee = -140:0.5:20; 
hip  = ones(length(knee),1)*mean_hip_flexion;
time = 0.01:0.01:length(knee)/100;
zero = zeros(length(knee),1);

if leg == 0;
    data_new = [time' zero zero zero zero zero zero hip zero zero knee' zero zero zero zero zero zero zero zero zero zero zero zero zero] ;
else
    data_new = [time' zero zero zero zero zero zero zero zero zero zero zero zero zero hip zero zero knee' zero zero zero zero zero zero] ;
end

% Write new mot file
q.data   = data_new;
q.labels = names;
fname    = 'C:\Users\u0125183\Box\PhD 1\Simulations Pendulum Test\Data\TD12/FakeMot_T2.mot';
write_motionFile(q, fname)

% %% Hip at zero 35 degrees
% addpath('C:\Users\u0125183\Documents\MATLAB\Sam/WriteMotFile')
% 
% % Read existing mot file
% [header, names, data, fpath] = SIMM_ReadMotion('C:\Users\u0125183\Documents\PhD 1\Pendulum test\Data\TD_05_07032019\Session 1\OpenSim\IK/Trial05.mot');
% 
% % new data 
% knee = -120:0.01:10; 
% hip  = ones(length(knee),1)*35;
% time = 0.01:0.01:length(knee)/100;
% zero = zeros(length(knee),1);
% 
% data_new = [time' zero zero zero zero zero zero zero zero zero zero zero zero zero hip zero zero knee' zero zero zero zero zero zero] 
% 
% % Write new mot file
% q.data   = data_new;
% q.labels = names;
% fname    = 'C:\Users\u0125183\Box\PhD 1\Simulations Pendulum Test\Spiermodel\OpenSim_data/FakeMot2.mot'
% write_motionFile(q, fname)
% 
% %% Hip at zero 90 degrees
% addpath('C:\Users\u0125183\Documents\MATLAB\Sam/WriteMotFile')
% 
% % Read existing mot file
% [header, names, data, fpath] = SIMM_ReadMotion('C:\Users\u0125183\Documents\PhD 1\Pendulum test\Data\TD_05_07032019\Session 1\OpenSim\IK/Trial05.mot');
% 
% % new data 
% knee = -120:0.01:10; 
% hip  = ones(length(knee),1)*90;
% time = 0.01:0.01:length(knee)/100;
% zero = zeros(length(knee),1);
% 
% data_new = [time' zero zero zero zero zero zero zero zero zero zero zero zero zero hip zero zero knee' zero zero zero zero zero zero] 
% 
% % Write new mot file
% q.data   = data_new;
% q.labels = names;
% fname    = 'C:\Users\u0125183\Box\PhD 1\Simulations Pendulum Test\Spiermodel\OpenSim_data/FakeMot3.mot'
% write_motionFile(q, fname)
% 
% %% Hip at zero 23 degrees (TD5)
% addpath('C:\Users\u0125183\Documents\MATLAB\Sam/WriteMotFile')
% 
% % Read existing mot file
% [header, names, data, fpath] = SIMM_ReadMotion('C:\Users\u0125183\Documents\PhD 1\Pendulum test\Data\TD_05_07032019\Session 1\OpenSim\IK/Trial05.mot');
% 
% % new data 
% knee = -120:0.01:10; 
% hip  = ones(length(knee),1)*90;
% time = 0.01:0.01:length(knee)/100;
% zero = zeros(length(knee),1);
% 
% data_new = [time' zero zero zero zero zero zero zero zero zero zero zero zero zero hip zero zero knee' zero zero zero zero zero zero] 
% 
% % Write new mot file
% q.data   = data_new;
% q.labels = names;
% fname    = 'C:\Users\u0125183\Box\PhD 1\Simulations Pendulum Test\Spiermodel\OpenSim_data/FakeMot_TD5.mot'
% write_motionFile(q, fname)
