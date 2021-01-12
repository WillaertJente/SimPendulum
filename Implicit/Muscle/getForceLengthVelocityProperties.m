function [Fpe, FMltilda, FMvtilda] = getForceLengthVelocityProperties(lMtilda, params, vMtilda, vMtildamax, kFpe)
%calculate Fpe, FMltilda
%   Fpe      = Normalized passive muscle force 
%   FMltilda = Normalized force-length multiplier
%   FMvtilda = Normalized force-velocity multiplier

% Fpe
e0  = 0.6;   kpe = 4;
t5  = exp(kpe * (lMtilda - kFpe*10) / e0);
Fpe = ((t5 - 0.10e1) - params.Fpparam(1)) / params.Fpparam(2);        %Fpe = musclepassiveforcelength(lMtilda);

% FMltilda
    % Active muscle force-length characteristics
b11 = params.Faparam(1); b21 = params.Faparam(2); b31 = params.Faparam(3); b41 = params.Faparam(4);
b12 = params.Faparam(5); b22 = params.Faparam(6); b32 = params.Faparam(7); b42 = params.Faparam(8);
b13 = 0.1;               b23 = 1;                 b33 = 0.5*sqrt(0.5);     b43 = 0;

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
    % active force-velocity
e1 = params.Fvparam(1); 
e2 = params.Fvparam(2);
e3 = params.Fvparam(3);
e4 = params.Fvparam(4); 

FMvtilda = e1*log((e2*vMtilda./vMtildamax+e3)+sqrt((e2*vMtilda./vMtildamax+e3).^2+1))+e4;

end

