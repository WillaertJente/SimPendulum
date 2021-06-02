function [params] = ImportParameters_forward(name)
% Import parameters for pendulum test
pathmain = pwd;
[pathTemp,~,~] = fileparts(pathmain);
[pathRepo,~,~] = fileparts(pathTemp);
info     = load([pathRepo '/Simpendulum\BK\' ,name]);

% Define parameters
params.m_tot   = info.(char(name)).m;               % Totall mass of subject
params.age     = info.(char(name)).age;             % Age of subject
params.z       = info.(char(name)).z;               % 18 if left leg, 11 if right leg
params.Nmr     = info.(char(name)).Nmr;             % First trial with MR
params.g       = 9.81;                              % Gravitational constant
params.length_tibia = info.(char(name)).lTibia; 
params.B       = info.(char(name)).B;               % Personalized damping
end