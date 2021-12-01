function [error] = CalculateMusculoSkeletalDynamics(q,qd,qdd,lMtilda,lM_projected,kFpe, vMtilda,a, data_exp, coeff_LMT_ma, params_OS, shift)

% Calculate Muscle tendon lengths and moment arms 
[lMT, MA] = CalculateMuscleTendonLengthAndMomentArms(q, data_exp, coeff_LMT_ma); 

% Calculate Tendon Force, muscle length and tendon length
[FT, lM, lT, fse, w ] = CalculateTendonForce(lMtilda,lM_projected, params_OS, lMT, shift);

% Get ForceLengthVelocity Relationships
[Fpe, FMltilda, FMvtilda] = getForceLengthVelocityRelation(lMtilda,kFpe, params_OS, vMtilda);

% FMce 
Fce = a.* FMltilda.* FMvtilda;      % FMce = fse.* lM ./(lMT-lT) - Fpe;

% Force equilibrium 
[error_force] = ForceEquilibrium(Fce, Fpe, lMT, lT, lM, fse); 

% Error lm and lM projected
error_musclegeometry = lM.^2 - w.^2 - lM_projected.^2; 

%Implicit skelet dynamics 
I = params_OS.inert.I_OS; 
m = params_OS.inert.mass_OS; 
l = params_OS.inert.lc_OS; 
g = 9.81; 
ma_ext = MA(1,:); 

error_skelet = qdd*I + m*g*l*cos(q) - FT.*ma_ext;  

% Total error
error = [error_force error_musclegeometry error_skelet];
end

