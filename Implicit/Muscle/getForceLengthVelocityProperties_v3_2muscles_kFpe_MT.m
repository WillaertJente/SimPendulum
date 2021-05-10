function [Fpe_ext,Fpe_flex, FMltilda_ext, FMltilda_flex, FMvtilda_ext, FMvtilda_flex] = getForceLengthVelocityProperties_v3_2muscles_kFpe_MT(lMtilda_ext, lMtilda_flex, params, vMtilda_ext, vMtilda_flex, vMtildamax_ext, vMtildamax_flex, kFpe_ext, kFpe_flex)
%calculate Fpe, FMltilda
%   Fpe      = Normalized passive muscle force 
%   FMltilda = Normalized force-length multiplier
%   FMvtilda = Normalized force-velocity multiplier

% Fpe
e0  = 0.6;   kpe = 4;
t5_ext  = exp(kpe * (lMtilda_ext - kFpe_ext*10) / e0);
Fpe_ext = ((t5_ext - 0.10e1) - params.Fpparam(1)) / params.Fpparam(2);        %Fpe = musclepassiveforcelength(lMtilda);

t5_flex  = exp(kpe * (lMtilda_flex - kFpe_flex*10) / e0);
Fpe_flex = ((t5_flex - 0.10e1) - params.Fpparam(1)) / params.Fpparam(2);      %Fpe = musclepassiveforcelength(lMtilda);

% FMltilda
    % Active muscle force-length characteristics Extensor
b11 = params.Faparam(1); b21 = params.Faparam(2); b31 = params.Faparam(3); b41 = params.Faparam(4);
b12 = params.Faparam(5); b22 = params.Faparam(6); b32 = params.Faparam(7); b42 = params.Faparam(8);
b13 = 0.1;               b23 = 1;                 b33 = 0.5*sqrt(0.5);     b43 = 0;

num3_ext = lMtilda_ext-b23;
den3_ext = b33+b43*lMtilda_ext;
FMtilda3_ext = b13*exp(-0.5*num3_ext.^2./den3_ext.^2);

num1_ext = lMtilda_ext-b21;
den1_ext = b31+b41*lMtilda_ext;
FMtilda1_ext = b11*exp(-0.5*num1_ext.^2./den1_ext.^2);

num2_ext = lMtilda_ext-b22;
den2_ext = b32+b42*lMtilda_ext;
FMtilda2_ext = b12*exp(-0.5*num2_ext.^2./den2_ext.^2);

FMltilda_ext = FMtilda1_ext+FMtilda2_ext+FMtilda3_ext;

% FMltilda
    % Active muscle force-length characteristics Flexor
b11 = params.Faparam(1); b21 = params.Faparam(2); b31 = params.Faparam(3); b41 = params.Faparam(4);
b12 = params.Faparam(5); b22 = params.Faparam(6); b32 = params.Faparam(7); b42 = params.Faparam(8);
b13 = 0.1;               b23 = 1;                 b33 = 0.5*sqrt(0.5);     b43 = 0;

num3_flex = lMtilda_flex-b23;
den3_flex = b33+b43*lMtilda_flex;
FMtilda3_flex = b13*exp(-0.5*num3_flex.^2./den3_flex.^2);

num1_flex = lMtilda_flex-b21;
den1_flex = b31+b41*lMtilda_flex;
FMtilda1_flex = b11*exp(-0.5*num1_flex.^2./den1_flex.^2);

num2_flex = lMtilda_flex-b22;
den2_flex = b32+b42*lMtilda_flex;
FMtilda2_flex = b12*exp(-0.5*num2_flex.^2./den2_flex.^2);

FMltilda_flex = FMtilda1_flex+FMtilda2_flex+FMtilda3_flex;

% FMvtilda
    % active force-velocity
e1 = params.Fvparam(1); 
e2 = params.Fvparam(2);
e3 = params.Fvparam(3);
e4 = params.Fvparam(4); 

% FMvtilda_ext  = e1*log((e2*vMtilda_ext./vMtildamax_ext+e3)+sqrt((e2*vMtilda_ext./vMtildamax_ext+e3).^2+1))+e4;
% FMvtilda_flex = e1*log((e2*vMtilda_flex./vMtildamax_flex+e3)+sqrt((e2*vMtilda_flex./vMtildamax_flex+e3).^2+1))+e4;

FMvtilda_ext  = e1*log((e2*vMtilda_ext+e3)+sqrt((e2*vMtilda_ext+e3).^2+1))+e4; 
FMvtilda_flex = e1*log((e2*vMtilda_flex+e3)+sqrt((e2*vMtilda_flex+e3).^2+1))+e4;
end

