function [FT, ma, dlMdt, err, lM, lT, Fce, Fpe, FM] = CalculateTendonForceAndMomentArm(x, params, lMtilda, a, sn, shift, vMtilda, lM_projected, coeff_LMT_ma, offset, kFpe, Fsrs)
%Function to calculate tendon force and moment arms
%   1. Calculate moment arm 
%   moment arm is calculated as the derivative of the muscle tendon length
%   2. Calculate tendon force
%       Tendon force = Peak isometric muscle force * Fse
%       Fse = exponential function of normalized tendon length
%       Normalized tendon length = tendon length/ tendon slack length
%       Tendon slack length  = input parameter
%       Tendon length = Muscle tendon length - muscle fiber length (+
%       andere factor)
%   3. Calculate derivative of lMtilda

%% Stap 1 LMT en Ma
offset   = offset*pi/180; 
m_offset = mean(offset);

lMT = coeff_LMT_ma(1) + coeff_LMT_ma(2)*(x+m_offset) + coeff_LMT_ma(3)*(x+m_offset).^2 + coeff_LMT_ma(4)*(x+m_offset).^3;
ma  = -coeff_LMT_ma(2) + -coeff_LMT_ma(3)*(x+m_offset) + -coeff_LMT_ma(4)*(x+m_offset).^2;

%% Stap 2
% w 
lMo    = params.MTparams(2,:); 
alphao = params.MTparams(4,:); 
w      = lMo.* sin(alphao);

% Hill type muscle model: geometric relationships
% lM (muscle fiber length) 
lM     = lMtilda.* lMo;

% lMT
% lMT    = interp1(params.LMT_qknee, params.LMT, x);

% lT (tendon length)
lT       = lMT - lM_projected;
% lT     = lMT - sqrt(abs(lM.^2 - w.^2)); 

% lTs (Tendon slack length)
lTs    = params.MTparams(3,:); 

% lTtilda
lTtilda= lT./lTs; 

% Fse
fse         = (exp(35*(lTtilda - 0.995)))/5-0.25 + shift;
%   fse(fse<0) = 0; 

% FMo
FMo = params.MTparams(1,:);

% Compute tendon force
FT = FMo.* fse; 

%% Stap 3
% vMtildamax
vMtildamax = params.MTparams(5,:); 

% Get force length velocity parameters
[Fpe, FMltilda, FMvtilda] = getForceLengthVelocityProperties(lMtilda, params, vMtilda, vMtildamax, kFpe)

% FMce 
Fce = a.* FMltilda.* FMvtilda;
% FMce = fse.* lM ./(lMT-lT) - Fpe;

% Compute dlMdt
dlMdt = vMtilda.* vMtildamax./ lMo;

% Compute Muscle force
FM  = Fce + Fpe;
% FM  = Fce;

% Force equilibrium
cos_alpha = (lMT-lT)./lM;
err       = FM.*cos_alpha - fse; 
end

