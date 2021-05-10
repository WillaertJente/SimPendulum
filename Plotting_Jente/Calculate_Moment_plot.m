%% Script to calculate ID moment and plot it
% Jente Willaert - 20 april 2021
clear all; close all ; clc
%% Import data
subject    = 'CP8';
damp_coeff = 0.0771;
trials     = 1:20;

pathmain       = pwd;
[pathTemp,~,~] = fileparts(pathmain);
[pathRepo,~,~] = fileparts(pathTemp);
path           = [pathRepo '/SimPendulum\Implicit\Muscle\Experimental data\' subject '\'];

map      = 'C:\Users\u0125183\Box\PhD 1\Simulations Pendulum Test\';
map_out  = ['C:\Users\u0125183\Box\PhD 1\Simulations Pendulum Test\Spiermodel\Plots/', subject,'/ID3_T']; 

for i = 1:length(trials)
    FileName = [map, 'Results/Result_',subject,'_T',num2str(trials(i)),'_Fmv.mat'];
    if exist(FileName,'File')
        % Simulated data
        load(FileName)
        % Params used for simulation
        params = ImportParameters(subject);
        addpath('C:\Users\u0125183\Documents\MATLAB\SimPendulum\Implicit\Muscle/MuscleModel');
        
        % Opensim model
        import org.opensim.modeling.*
        model_path = [path,subject,'_ScaledModel_ScaledForces.osim'];                                 % if cp = CPModel_Scaled.osim
        osimModel  = Model(model_path);
        
        % Inertial parameters (tibia)
        % Massa
        bodies         = osimModel.getBodySet();
        if params.z == 18
            tibia      = bodies.get('tibia_l');
            talus      = bodies.get('talus_l');
            calcn      = bodies.get('calcn_l');
            toes       = bodies.get('toes_l');
        else
            tibia      = bodies.get('tibia_r');
            talus      = bodies.get('talus_r');
            calcn      = bodies.get('calcn_r');
            toes       = bodies.get('toes_r');
        end
        params.mass_OStibia = tibia.getMass();
        params.mass_OScalcn = calcn.getMass();
        params.mass_OS = params.mass_OStibia + params.mass_OScalcn;
        
        % COM
        com_tibia = ArrayDouble.createVec3(0);
        tibia.getMassCenter(com_tibia);
        com_calcn = ArrayDouble.createVec3(0);
        calcn.getMassCenter(com_calcn);
        com_foot_in_tibia = abs(com_calcn.get(1)) + params.length_tibia;
        params.lc_OS = (abs(com_tibia.get(1)) * tibia.getMass() + com_foot_in_tibia*calcn.getMass())/params.mass_OS;
        
        % Inertia
        inertia_tibia = Mat33(0);
        tibia.getInertia(inertia_tibia);
        inertia_calcn = Mat33(0);
        calcn.getInertia(inertia_calcn);
        params.I_OS = inertia_tibia.get(0,0) + inertia_calcn.get(0,0) + tibia.getMass() * com_tibia.get(1) ^ 2 + calcn.getMass() * com_foot_in_tibia ^ 2;
        
        % ID from openSim
        ID_OS   = importdata([map,'Data/',subject,'/T',num2str(trials(i)),'_inverse_dynamics.sto']);
        start_f = find(ID_OS.data(:,1)==tvect(1))
        if params.z == 18;
            col = find(strcmp(ID_OS.colheaders,'knee_angle_l_moment'));
        else
            col = find(strcmp(ID_OS.colheaders,'knee_angle_r_moment'));
        end
        
        %% Calculate ID moment
        xd = [diff(sol_x) 0];
        
        %         ID_msc_inert          = -params.mass_OS*params.g*params.lc_OS*cos(sol_x)+ sol_FT_ext.*sol_ma_ext + sol_FT_flex.*sol_ma_flex;
        %         ID_msc_inert_damp     = -params.mass_OS*params.g*params.lc_OS*cos(sol_x)+ sol_FT_ext.*sol_ma_ext + sol_FT_flex.*sol_ma_flex - damp_coeff*xd;
        %         ID_msc_inert_damp_act = -params.mass_OS*params.g*params.lc_OS*cos(sol_x)+ sol_FT_ext.*sol_ma_ext + sol_FT_flex.*sol_ma_flex + sol_act - damp_coeff*xd;
        
        %ID_msc_inert          =  sol_FT_ext.*sol_ma_ext + sol_FT_flex.*sol_ma_flex;
        %ID_msc_inert_damp     =  sol_FT_ext.*sol_ma_ext + sol_FT_flex.*sol_ma_flex - damp_coeff*xd;
        ID_msc_inert_damp_act  =  sol_FT_ext.*sol_ma_ext + sol_FT_flex.*sol_ma_flex + sol_act - damp_coeff*xd;
        %% Plot experiemntal, simulated and ID-moment data
        figure()
        % Experimental and simulated data
        subplot(211)
        plot(tvect,q_exp*180/pi,'k','LineWidth',1.5); hold on;
        plot(tvect,sol_x*180/pi,'r','LineWidth',1.5); hold on;
        ylabel('Knee angle{\circ})'); xlabel('Time (frames)')
        box off; legend('Exdxp','Sim')
        title([subject, ' T',num2str(trials(i))])
        
        % ID calculated by OS
        subplot(212)
        plot(ID_OS.data(start_f:end,1),ID_OS.data(start_f:end,col),'k','LineWidth',1.5); hold on
        plot(tvect,ID_msc_inert_damp_act,'r','LineWidth',1.5); hold on
        ylabel('Nm'); box off; legend({'OpenSim','Sim'});
        
%         % ID muscle + inertia
%         subplot(513)
%         plot(tvect,ID_msc_inert,'k','LineWidth',1.5);
%         ylabel('Nm'); legend('M+I'); box off;
%         
%         % ID muscle + inertia + damping
%         subplot(514)
%         plot(tvect,ID_msc_inert_damp,'k','LineWidth',1.5);
%         ylabel('Nm'); legend('M+I+D'); box off;
%         
%         % ID muscle + inertia + damping + actuator
%         subplot(515)
%         plot(tvect,ID_msc_inert_damp_act,'k','LineWidth',1.5);
%         ylabel('Nm'); legend('M+I+D+A'); box off ;
        
        saveas(gcf, [map_out,num2str(i),'.png'])
    end
end