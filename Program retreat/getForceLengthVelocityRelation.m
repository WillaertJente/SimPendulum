function [Fpe, FMltilda, FMvtilda] = getForceLengthVelocityRelation(lMtilda,kFpe, params_OS, vMtilda)
%calculate Fpe, FMltilda
%   Fpe      = Normalized passive muscle force 
%   FMltilda = Normalized force-length multiplier
%   FMvtilda = Normalized force-velocity multiplier

% Fpe
e0      = 0.6;   kpe = 4;
t5      = exp(kpe * (lMtilda - kFpe*10) / e0);          % Change it if we have a flexor kFpe as well
Fpe_ext = ((t5 - 0.10e1) - params_OS.Fp(1)) / params_OS.Fp(2);        %Fpe = musclepassiveforcelength(lMtilda);
Fpe     = Fpe_ext; 

% FMltilda
% Active muscle force-length characteristics
b11 = params_OS.Fa(1); b21 = params_OS.Fa(2); b31 = params_OS.Fa(3); b41 = params_OS.Fa(4);
b12 = params_OS.Fa(5); b22 = params_OS.Fa(6); b32 = params_OS.Fa(7); b42 = params_OS.Fa(8);
b13 = 0.1;             b23 = 1;               b33 = 0.5*sqrt(0.5);   b43 = 0;

num3 = lMtilda-b23;
den3 = b33+b43*lMtilda;
FMtilda3 = b13*exp(-0.5*num3.^2./den3.^2);

num1 = lMtilda-b21;
den1 = b31+b41*lMtilda;
FMtilda1 = b11*exp(-0.5*num1.^2./den1.^2);

num2 = lMtilda-b22;
den2 = b32+b42*lMtilda;
FMtilda2 = b12*exp(-0.5*num2.^2./den2.^2);

FMltilda = FMtilda1+FMtilda2+FMtilda3;

% FMvtilda
vMtildamax = params_OS.MT(5,:); 
% active force-velocity
e1 = params_OS.Fv(1); 
e2 = params_OS.Fv(2);
e3 = params_OS.Fv(3);
e4 = params_OS.Fv(4); 

FMvtilda_ext  = e1*log((e2*vMtilda./vMtildamax(1)+e3)+sqrt((e2*vMtilda./vMtildamax(1)+e3).^2+1))+e4; % extensor
FMvtilda_flex = e1*log((e2*vMtilda./vMtildamax(2)+e3)+sqrt((e2*vMtilda./vMtildamax(2)+e3).^2+1))+e4; % flexor, change vMtilda to flexor 
FMvtilda      = FMvtilda_ext; 
end

