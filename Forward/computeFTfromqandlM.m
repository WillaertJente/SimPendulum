function FT = computeFTfromqandlM(MTparams,coeff,q,lMtilde)
% Friedl De Groote
% June 8, 2021

FMo    = MTparams(1,:);
lMo    = MTparams(2,:);
lTs    = MTparams(3,:);
alphao = MTparams(4,:);
m_offset = 0; % @Jente - you have this variable in way too many places
lMT = computeLengthMomentArm(coeff,q+m_offset);

lM = lMtilde.*lMo;
w  = lMo.*sin(alphao);
lT = lMT - sqrt(abs(lM.^2 - w.^2));
lTtilda = lT./lTs;

% Compute tendon force
fse = (exp(35*(lTtilda - 0.995)))/5-0.25;
fse(fse<0) = 0;
FT = FMo.*fse;
end

