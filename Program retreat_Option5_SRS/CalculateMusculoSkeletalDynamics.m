function [error] = CalculateMusculoSkeletalDynamics(q,qd,qdd,lMtilda, lM_projected, kFpe, vMtilda, a_ext, a_flex, data_exp, coeff_LMT_ma, params_OS, shift, B, info, Fsrs, dFsrsdt)

% Calculate Muscle tendon lengths and moment arms 
[lMT, MA] = CalculateMuscleTendonLengthAndMomentArms(q, data_exp, coeff_LMT_ma); 

% Calculate Tendon Force, muscle length and tendon length
[FT, lM, lT, fse, w ] = CalculateTendonForce(lMtilda, lM_projected, params_OS, lMT, shift);

% Get ForceLengthVelocity Relationships
[Fpe, FMltilda, FMvtilda] = getForceLengthVelocityRelation(lMtilda, kFpe, params_OS, vMtilda);

% Fsrs 

% % Calculate Fsrs 
% [Fsrs, dFsrs2dt_cal] = CalculateSRS(info, data_exp, lMtilda, a_ext, Fsrs2, FMltilda); 

% % Fsrs_error
% N_1 = data_exp.N_1; 
% error_Fsrs = dFsrs2dt_cal - dFsrs2dt; 

% FMce 
Fce_ext  = a_ext.* FMltilda(1,:).* FMvtilda(1,:) + Fsrs;      % FMce = fse.* lM ./(lMT-lT) - Fpe;
Fce_flex = a_flex.* FMltilda(2,:).* FMvtilda(2,:); 
Fce      = [Fce_ext; Fce_flex];

% Force equilibrium 
[error_force] = ForceEquilibrium(Fce, Fpe, lMT, lT, lM, fse); 

% Error lm and lM projected
error_musclegeometry_ext  = lM(1,:).^2 - w(1,:).^2 - lM_projected(1,:).^2; 
error_musclegeometry_flex = lM(2,:).^2 - w(2,:).^2 - lM_projected(2,:).^2;
error_musclegeometry      = [error_musclegeometry_ext error_musclegeometry_flex]; 

%Implicit skelet dynamics 
I = params_OS.inert.I_OS; 
m = params_OS.inert.mass_OS; 
l = params_OS.inert.lc_OS; 
g = 9.81; 
ma_ext = MA(1,:); 
ma_flex = MA(2,:); 

error_skelet = qdd*I + m*g*l*cos(q) - FT(1,:).*ma_ext - FT(2,:).*ma_flex + B*qd ;  

% Total error
error = [error_force error_musclegeometry error_skelet error_Fsrs];
end

