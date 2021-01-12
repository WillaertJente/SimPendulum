%% calculate tendon force and moment arm voor bepaalde q_exp

%% Input
% Experimentele data
s.nu   = 'TD5';                                                                               % subject number/ name
s.tr   = 1;                                                                                 % subject trials (number of trials)
path   = 'C:\Users\u0125183\Box\PhD 1\Pendulum test\Data\TD_05_07032019\Session 1\OpenSim';   % Path to opensim model (scaled)
params = ImportParameters(s.nu);
j      = 1;

[q_exp_r, t_exp, t_span, on_srs] = LoadExpData(s.nu,s.tr,j,0, params);      % 1 if you want to plot experimental data
[q_exp, qdot_exp,N, tvect, dt]   = ExpDatAtDiscrTime(t_span,t_exp,q_exp_r);

x = q_exp;

% Params
addpath('MuscleModel');
import org.opensim.modeling.*
model_path = [path,'/TDModel_Scaled_MusclAnal.osim'];                                 % if cp = CPModel_Scaled.osim
osimModel  = Model(model_path);
[params_Muscle,lOpt,L_TendonSlack,Fiso,PennationAngle]=ReadMuscleParameters(model_path,{'rect_fem_r'});
params_Muscle(2,1) = params_Muscle(2,1);
params_Muscle(3,1) = params_Muscle(3,1);
params.MTparams    = params_Muscle;

load Fvparam.mat
Fvparam = [-0.2158 -32.5966 -1.1241 0.9126];
Fvparam(1)       = 1.475*Fvparam(1); Fvparam(2) = 0.25*Fvparam(2);
Fvparam(3)       = Fvparam(3) + 0.75; Fvparam(4) = Fvparam(4) - 0.027;
params.Fvparam   = Fvparam;


load Faparam.mat
params.Faparam   = Faparam;


e0  = 0.6; kpe = 4; t50  = exp(kpe * (0.2 - 0.10e1) / e0);
pp1 = (t50 - 0.10e1); t7 = exp(kpe); pp2 = (t7 - 0.10e1);
params.Fpparam = [pp1;pp2];

% lMtilda
lMtilda = 0.9;

% a
a = 0.1;

% shift
shift = 0;

% vMtilda
vMtilda = -0.3;

% lMprojected
lMo      = params.MTparams(2,:)';
alphao   = params.MTparams(4,:)';
w        = lMo.*sin(alphao);            % in radians
lM       = lMtilda.*lMo;
lM_projected   = sqrt(lM.^2 - w^2);


% Coeff LMT_ma
map_MA         = 'C:\Users\u0125183\Box\PhD 1\Simulations Pendulum Test\Spiermodel\OpenSim_data\Muscle Analysis/';
[coeff_LMT_ma] = DefineLMTCoefficients(map_MA, s.nu);

%% Calculate tendon force and moment arm
[FT, ma, dlMdt, err, lM] = CalculateTendonForceAndMomentArm(x, params, lMtilda, a, s.nu, shift, vMtilda, lM_projected, coeff_LMT_ma)

lMT      = coeff_LMT_ma(1) + coeff_LMT_ma(2)*q_exp + coeff_LMT_ma(3)*q_exp.^2 + coeff_LMT_ma(4)*q_exp.^3
lT       = lMT - lM_projected;
lTs      = params.MTparams(3,:);
