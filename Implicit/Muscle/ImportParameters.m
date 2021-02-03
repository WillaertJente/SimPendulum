function [params] = ImportParameters(name)
% Import parameters for pendulum test
pathmain = pwd;
[pathTemp,~,~] = fileparts(pathmain);
[pathRepo,~,~] = fileparts(pathTemp);
info     = load([pathRepo '\BK\' ,name]);

% Define parameters
params.m_tot   = info.(char(name)).m;               % Totall mass of subject
params.age     = info.(char(name)).age;             % Age of subject
%params.lc      = info.(char(name)).lc  ;            % Distance knee - COM tibia
%params.l       = info.(char(name)).l;               % Leg length
params.z       = info.(char(name)).z;               % 18 if left leg, 11 if right leg
params.Nmr     = info.(char(name)).Nmr;             % First trial with MR
%m_ratio        = 0.00137 * params.age + (0.03809+0.0187);   % ratio to calculate leg mass
%params.m       = m_ratio * params.m_tot;            % Mass of the leg
params.g       = 9.81;                              % Gravitational constant
%params.RG      = params.l * 0.416;                  % Radius of gyration (Winter 2009)
%params.SE      = (0.0055 + 0.0028)*params.m_tot;    % SE used in equation to estimate leg mass
params.length_tibia = info.(char(name)).lTibia; 

end



% m_ratio_voet  = (0.00015*age + 0.0187);
% m_ratio_shank = (0.00122*age + 0.03809);
% m_SE_voet     = 0.0028;
% m_SE_shank    = 0.0055; 
% RG_voet       = -0.00203*age + 0.5022;
% RG_shank      = -0.00224*age + 0.5307;
% RG_SE_voet    = 0.0161;
% RG_SE_shank   = 0.0112; 
% a             = -0.00186*age + 0.4351;
% a_SE          = 0.0123;
% b             = l;
% theta         = 44*pi/180;
% lc_shank      = -0.003*age + 0.4526;
% lc_SE_shank   = 0.0123; 
% d             = sqrt(a*a + b*b -2*a*b*cos(theta));
% 
% m       = (m_ratio_voet + m_ratio_shank)* m_tot;
% g       = 9.81;                                      % gravitational constant
%            
% m_voet  = m_ratio_voet*m_tot;
% m_ob    = m_ratio_shank*m_tot;
% I       = (m_voet*RG_voet*RG_voet) + (m_voet*d*d) + (m_ob*RG_shank*RG_shank)+ (m_ob * lc_shank * lc_shank)  ;
% 
% disp(I)