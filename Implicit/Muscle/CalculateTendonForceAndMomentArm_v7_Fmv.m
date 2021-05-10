function [FT_ext,FT_flex, ma_ext, ma_flex, dlMdt_ext, dlMdt_flex, err_ext, err_flex, lM_ext, lM_flex, lT_ext, lT_flex, Fce_ext, Fce_flex, Fpe_ext, Fpe_flex, FM_ext, FM_flex, Fsrs, FMltilda_ext, FMltilda_flex, FMvtilda_ext] = CalculateTendonForceAndMomentArm_v7_Fmv(x, params, lMtilda_ext, lMtilda_flex, a_ext_0,a_flex, shift, vMtilda_ext,vMtilda_flex, lM_projected_ext, lM_projected_flex, coeff_LMT_ma_ext, coeff_LMT_ma_flex, m_offset, kFpe_ext, kFpe_flex, N_1, N, Fsrs_d, Fsrs2, a_ext)
%Function to calculate tendon force and moment arms with SRS 
%   1. Calculate moment arm
%   ma is derivative of LMT
%   2. Calculate tendon force
%       Tendon force = Peak isometric muscle force * Fse
%       Fse = exponential function of normalized tendon length
%       Normalized tendon length = tendon length/ tendon slack length
%       Tendon slack length  = input parameter
%       Tendon length = Muscle tendon length - muscle fiber length (+
%       andere factor)
%   3. Calculate derivative of lMtilda

%% Stap 1 LMT en Ma
lMT_ext = coeff_LMT_ma_ext(1) + coeff_LMT_ma_ext(2)*(x+m_offset) + coeff_LMT_ma_ext(3)*(x+m_offset).^2 + coeff_LMT_ma_ext(4)*(x+m_offset).^3 + coeff_LMT_ma_ext(5)*(x+m_offset).^4 + coeff_LMT_ma_ext(6)*(x+m_offset).^5;
ma_ext  = -coeff_LMT_ma_ext(2) + -coeff_LMT_ma_ext(3)*(x+m_offset) + -coeff_LMT_ma_ext(4)*(x+m_offset).^2 + -coeff_LMT_ma_ext(5)*(x+m_offset).^3 + -coeff_LMT_ma_ext(6)*(x+m_offset).^4;

lMT_flex = coeff_LMT_ma_flex(1) + coeff_LMT_ma_flex(2)*(x+m_offset) + coeff_LMT_ma_flex(3)*(x+m_offset).^2 + coeff_LMT_ma_flex(4)*(x+m_offset).^3 + coeff_LMT_ma_flex(5)*(x+m_offset).^4 + coeff_LMT_ma_flex(6)*(x+m_offset).^5;
ma_flex  = -coeff_LMT_ma_flex(2) + -coeff_LMT_ma_flex(3)*(x+m_offset) + -coeff_LMT_ma_flex(4)*(x+m_offset).^2 + -coeff_LMT_ma_flex(5)*(x+m_offset).^3 + -coeff_LMT_ma_flex(6)*(x+m_offset).^4;

%% Stap 2
lMo_ext    = params.MTparams_ext(2,:);
lMo_flex   = params.MTparams_flex(2,:);

% Hill type muscle model: geometric relationships
% lM (muscle fiber length)
lM_ext     = lMtilda_ext.* lMo_ext;
lM_flex    = lMtilda_flex.* lMo_flex;

% lT (tendon length)
lT_ext     = lMT_ext  - lM_projected_ext;
lT_flex    = lMT_flex - lM_projected_flex;

% lTs (Tendon slack length)
lTs_ext    = params.MTparams_ext(3,:);
lTs_flex   = params.MTparams_flex(3,:);

% lTtilda
lTtilda_ext  = lT_ext./lTs_ext;
lTtilda_flex = lT_flex./lTs_flex;

% Fse
fse_ext      = (exp(35*(lTtilda_ext - 0.995)))/5-0.25 + shift;
fse_flex     = (exp(35*(lTtilda_flex - 0.995)))/5-0.25 + shift;

% FMo
FMo_ext  = params.MTparams_ext(1,:);
FMo_flex = params.MTparams_flex(1,:);

% Compute tendon force
FT_ext  = FMo_ext.* fse_ext;
FT_flex = FMo_flex.* fse_flex;

%% Stap 3
% vMtildamax
vMtildamax_ext  = params.MTparams_ext(5,:);
vMtildamax_flex = params.MTparams_flex(5,:);

% Get force length velocity parameters
% Friedl - would be better to use vector notation and not create new
% functions every time you change number of muscles.
[Fpe_ext,Fpe_flex, FMltilda_ext, FMltilda_flex, FMvtilda_ext, FMvtilda_flex] = getForceLengthVelocityProperties_v3_2muscles_kFpe(lMtilda_ext, lMtilda_flex, params, vMtilda_ext, vMtilda_flex, vMtildamax_ext, vMtildamax_flex, kFpe_ext, kFpe_flex);

FMvtilda_ext(1:N_1) = 1; 

% Tsrs
kSRS = 280;
dLm  = lMtilda_ext - lMtilda_ext(1);        % Stretch 

Fsrs     =(0.5*tanh(1000*(-dLm(1:N_1)+5.7*10^(-3)))+0.5).*dLm(1:N_1).*FMltilda_ext(1:N_1)*a_ext_0*kSRS + ...
    (0.5*tanh(1000*(dLm(1:N_1) - 5.7*10^(-3)))+0.5)*5.7*10^(-3).*a_ext_0.*FMltilda_ext(1:N_1)*kSRS;    

% FMce
Fce_ext  = a_ext.* FMltilda_ext.* FMvtilda_ext + [Fsrs Fsrs2]; 
Fce_flex = a_flex.* FMltilda_flex.* FMvtilda_flex; 

% Compute dlMdt
% dlMdt_ext  = vMtilda_ext.* vMtildamax_ext./ lMo_ext;
% dlMdt_flex = vMtilda_flex.* vMtildamax_flex./ lMo_flex;
dlMdt_ext  = vMtilda_ext.* vMtildamax_ext;
dlMdt_flex = vMtilda_flex.* vMtildamax_flex;

% Compute Muscle force
% Friedl - do not allow the passive forces to become larger than 1.5
Fpe_ext_lim  = (- 0.5*tanh(10*(Fpe_ext-1.5)) + 0.5).*Fpe_ext + (0.5*tanh(10*(Fpe_ext-1.5)) + 0.5)*1.5;
Fpe_flex_lim = (- 0.5*tanh(10*(Fpe_flex-1.5)) + 0.5).*Fpe_flex + (0.5*tanh(10*(Fpe_flex-1.5)) + 0.5)*1.5;
FM_ext   = Fce_ext + Fpe_ext_lim;
FM_flex  = Fce_flex + Fpe_flex_lim;

% FM_ext   = Fce_ext + Fpe_ext;
% FM_flex  = Fce_flex + Fpe_flex;

% Force equilibrium
cos_alpha_ext = (lMT_ext-lT_ext)./lM_ext;
err_ext       = FM_ext.*cos_alpha_ext - fse_ext;

cos_alpha_flex = (lMT_flex-lT_flex)./lM_flex;
err_flex       = FM_flex.*cos_alpha_flex - fse_flex;
end

