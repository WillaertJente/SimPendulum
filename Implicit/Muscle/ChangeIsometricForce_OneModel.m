%%% Change isometric force of the subjects
%=========================================
%Sam Van Rossom - 07-09-2017

clear all, close all, clc

%% Input 
%-------

path_gen = 'C:\Users\u0125183\Box\PhD 1\Simulations Pendulum Test\Data\CP8\CP8_ScaledModel_Torso.osim';
% InfoFile = 'J:\GBW-0301_HumanMovementBiomechanics\SimCP\OpenSimWorkflow_TestSam\SubjectsSam2.xlsx';
% root = 'J:\GBW-0301_HumanMovementBiomechanics\SimCP\Model_Database\BOTOX'; 
% rootxml = 'J:\GBW-0301_HumanMovementBiomechanics\SimCP\OpenSimWorkflow_TestSam\ReferenceSetup'; % root directory where your reference OpenSim setup files are
    
%% Scale model
%-------------

%Read info from Excel
%--------------------
% [num,txt,raw] = xlsread(InfoFile,'Static');

% load the API
%-------------
import org.opensim.modeling.*

%load generic model
%------------------
GenModel = Model(path_gen); 
GenModel.initSystem;
GenMuscles = GenModel.getMuscles(); 

% Count the muscles
%-----------------
nMuscles = GenMuscles.getSize();


% for p = 2:size(raw,1)
 mass =33.8; 
   scaleFactor =((mass/75.16).^(2/3));
   
   model_in = 'C:\Users\u0125183\Box\PhD 1\Simulations Pendulum Test\Data\CP8\CP8_ScaledModel_ScaledForces_Torso.osim'; %fullfile(root,raw{p,1},'static',[raw{p,1} '.osim']); 
   
   MyModel = Model(model_in); 
   MyModel.initSystem; 
   
   MyMuscles = MyModel.getMuscles(); 
%    nMuscles = MyMuscles.getsize(); 
   
   for m = 0:nMuscles - 1; 
      GenMuscle = GenMuscles.get(m); 
      ori_force =  GenMuscle.get_max_isometric_force(); 
      scaled_force = ori_force*scaleFactor;
      
      MyMuscle = MyMuscles.get(m); 
      MyMuscle.set_max_isometric_force(scaled_force);
   end 
   MyModel.print(model_in);  
   
   
   %Create actuators file
%    modelFile = [root '\' raw{p,1} '\static\' raw{p,1} '.osim'];
%    
%    ForceSetFile_org = [rootxml '\SO_Actuators.xml'];
%    ForceSetFile = [root '\' raw{p,1} '\static\' raw{p,1} '_SO_Actuators.xml'];
%    copyfile(ForceSetFile_org,ForceSetFile,'f');
%    
%    CreateActuatorsFile(ForceSetFile_org,ForceSetFile,modelFile);
             
% end



