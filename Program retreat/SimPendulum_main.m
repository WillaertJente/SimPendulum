%% SimPendulum_main

%% Input 
% Change here subject and trial
info.subj   = 'TD5';           % Subject name                                                                    
info.trial  = 2;               % Trial number                                                                     
info.option = '';              % Name to save results 

%% Import parameters and experimental data
 
% Path info -  Path to model and experimental data 
pathmain = pwd;
[pathTemp,~,~] = fileparts(pathmain);
[pathRepo,~,~] = fileparts(pathTemp);
info.path      = [pathTemp '\Implicit\Muscle\Experimental data\' info.subj '\'];        

% Import subject parameters 
params_subject = ImportSubjectParameters(info);    % Input parameters (mtot,lc, l, age, m, RG, SE, Nmr, z)

% Import experimental data 
bool_plot = 1;
dt_spline = 0.005; 
data_exp  = ImportExperimentalData(info, bool_plot, params_subject, dt_spline);
 
% Define phases of pendulum (initial state, end of first swing) 
bool_plot = 1;
[data_exp.x0, data_exp.N_1] = PendulumPhases(data_exp, bool_plot);  
         
%% Import OpenSim Model 

