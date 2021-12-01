function [output] = RecalculateOutcomes(R)
% Function to recalculate Forces and Lengths

% Input
x        = R.x; 
a        = R.a; 
lMtilda  = R.lMtilda; 
kFpe     = R.kFpe; 
vMtilda  = R.vMtilda; 
offset   = R.exp.offset; 
m_offset = mean(offset);
lTs      = R.OS.MT(3,:); 
lMo      = R.OS.MT(2,:); 
FMo      = R.OS.MT(1,:);
shift    = R.shift; 
coeff_LMT_ma = R.coeff; 
lM_projected = R.lMprojected;  
vMtildamax = R.OS.MT(5,:); 

% Muscle Tendon Lengths
lMT_ext     = coeff_LMT_ma(1,1) + coeff_LMT_ma(2,1)*(x+m_offset) + coeff_LMT_ma(3,1)*(x+m_offset).^2 + coeff_LMT_ma(4,1)*(x+m_offset).^3;
output.lMT_ext = lMT_ext; 

% Moment Arms
MA_ext  = -coeff_LMT_ma(2,1) + -coeff_LMT_ma(3,1)*(x+m_offset) + -coeff_LMT_ma(4,1)*(x+m_offset).^2;
output.MA_ext = MA_ext; 

% lT 
lT_ext   = lMT_ext - lM_projected; 
output.lT_ext = lT_ext; 

% lTtilda
lTtilda_ext  = lT_ext./lTs(1); 
output.lTtilda_ext = lTtilda_ext; 

% lM 
lM_ext  = lMtilda.* lMo(1);
output.lM_ext = lM_ext; 

% Fse
fse_ext      = (exp(35*(lTtilda_ext - 0.995)))/5-0.25 + shift;
output.fse_ext = fse_ext; 

% FT
FT_ext       = FMo(1).* fse_ext; 
output.FT_ext = FT_ext; 

% Fpe
e0      = 0.6;   kpe = 4;
t5      = exp(kpe * (lMtilda - kFpe*10) / e0);          % Change it if we have a flexor kFpe as well
Fpe_ext = ((t5 - 0.10e1) - R.OS.Fp(1)) / R.OS.Fp(2);        %Fpe = musclepassiveforcelength(lMtilda);
Fpe     = Fpe_ext; 
output.Fpe_ext = Fpe_ext; 

% FMltilda
% Active muscle force-length characteristics
b11 = R.OS.Fa(1); b21 = R.OS.Fa(2); b31 = R.OS.Fa(3);       b41 = R.OS.Fa(4);
b12 = R.OS.Fa(5); b22 = R.OS.Fa(6); b32 = R.OS.Fa(7);       b42 = R.OS.Fa(8);
b13 = 0.1;        b23 = 1;          b33 = 0.5*sqrt(0.5);    b43 = 0;

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
output.FMltilda = FMltilda;

% FMvtilda
% active force-velocity
e1 = params_OS.Fv(1); 
e2 = params_OS.Fv(2);
e3 = params_OS.Fv(3);
e4 = params_OS.Fv(4); 

FMvtilda_ext  = e1*log((e2*vMtilda./vMtildamax(1)+e3)+sqrt((e2*vMtilda./vMtildamax(1)+e3).^2+1))+e4; % extensor
FMvtilda_flex = e1*log((e2*vMtilda./vMtildamax(2)+e3)+sqrt((e2*vMtilda./vMtildamax(2)+e3).^2+1))+e4; % flexor, change vMtilda to flexor 
FMvtilda      = FMvtilda_ext; 
output.FMvtilda = FMvtilda; 

% Fce
Fce = a.* FMltilda.* FMvtilda; 
output.Fce = Fce; 

end

