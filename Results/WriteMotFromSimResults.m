addpath('C:\Users\u0125183\Documents\MATLAB\NeuromechanicsToolkit-main\ReadWrite/')

sim = load('TD5_T1_Opt12_IG_set5_Bounds_Set1.mat');

colnames = {'time'	'pelvis_tilt'	'pelvis_list'	'pelvis_rotation'...
    'pelvis_tx'	'pelvis_ty'	'pelvis_tz'	'hip_flexion_r'	'hip_adduction_r'...
    'hip_rotation_r'	'knee_angle_r'	'ankle_angle_r'	'subtalar_angle_r'...
    'mtp_angle_r'	'hip_flexion_l'	'hip_adduction_l'	'hip_rotation_l'...
    'knee_angle_l'	'ankle_angle_l'	'subtalar_angle_l'	'mtp_angle_l'...
    'lumbar_extension'	'lumbar_bending'	'lumbar_rotation'}; 

data = zeros(length(R.exp.tspline),length(colnames));
data(:,1)  = R.exp.tspline;
data(:,18) = R.x;
dataMatrix = data;

filename= ('C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Implicit\Muscle\Experimental data\New/CP10/T13_FakeMotFromBK.mot'); 

generateMotFile(dataMatrix, colnames, filename)