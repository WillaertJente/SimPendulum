function [dlMdt,FT] = FiberLengthOde_flex(a,lMtilda,lMT,MTparams, Fvparam, Fpparam, Faparam, kFpe_flex)
% Try to vectorize it
FMo = MTparams(1,:);              %Peak isometric muscle force
lMo = MTparams(2,:);              %Optimal muscle fiber length
lTs = MTparams(3,:);              %Tendon slack length 
alphao = MTparams(4,:);           %Pennation angle at optimal fiber length
vMmax = MTparams(5,:);            %Maximal fiber velocity 

lM = lMtilda.*lMo;
w = lMo.*sin(alphao);
lT = lMT - sqrt(abs(lM.^2 - w.^2));
lTtilda = lT./lTs;

% Compute tendon force
fse = (exp(35*(lTtilda - 0.995)))/5-0.25;
fse(fse<0) = 0;
FT = FMo.*fse;

% Gaussians
b11 = Faparam(1);
b21 = Faparam(2);
b31 = Faparam(3);
b41 = Faparam(4);
b12 = Faparam(5);
b22 = Faparam(6);
b32 = Faparam(7);
b42 = Faparam(8);

b13 = 0.1;
b23 = 1;
b33 = 0.5*sqrt(0.5);
b43 = 0;
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

% Calculate passive force length
% Fpe = musclepassiveforcelength(lMtilda);
e0 = 0.6;
kpe = 4;
t5 = exp(kpe * (lMtilda - kFpe_flex *10) / e0);
Fpe = ((t5 - 0.10e1) - Fpparam(1)) / Fpparam(2);

FMce = fse.* lM ./(lMT-lT) - Fpe;
FMvtilda = FMce./(a.*FMltilda);
FMvtilda(FMvtilda<0) = 0;
FMvtilda(FMvtilda>1.8) = 1.8;

% FMvtilda = muscleforcevelocity(vMtilda);
% load Fvparam
e1 = Fvparam(1);
e2 = Fvparam(2);
e3 = Fvparam(3);
e4 = Fvparam(4);

vMtilda = 1/e2*(sinh((FMvtilda-e4)/e1)-e3);
dlMdt = vMtilda .* vMmax ./ lMo;

return