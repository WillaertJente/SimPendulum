function [Fsrs] = CalculateFsrsF1Forward(lMtilda_forward_F1,lMtilda_init, a_ext, kSRS, params_OS)
%Function to calculate final value of SRS at end of F1 
%   Needed as input for the forward calculation of F2 

% Stretch at end of F1
stretch = lMtilda_forward_F1-lMtilda_init; 

% FMltilda
% Active muscle force-length characteristics
b11 = params_OS.Fa(1); b21 = params_OS.Fa(2); b31 = params_OS.Fa(3); b41 = params_OS.Fa(4);
b12 = params_OS.Fa(5); b22 = params_OS.Fa(6); b32 = params_OS.Fa(7); b42 = params_OS.Fa(8);
b13 = 0.1;             b23 = 1;               b33 = 0.5*sqrt(0.5);   b43 = 0;

% ext 
num3_ext     = lMtilda_forward_F1-b23;
den3_ext     = b33+b43*lMtilda_forward_F1;
FMtilda3_ext = b13*exp(-0.5*num3_ext.^2./den3_ext.^2);

num1_ext     = lMtilda_forward_F1-b21;
den1_ext     = b31+b41*lMtilda_forward_F1;
FMtilda1_ext = b11*exp(-0.5*num1_ext.^2./den1_ext.^2);

num2_ext     = lMtilda_forward_F1-b22;
den2_ext     = b32+b42*lMtilda_forward_F1;
FMtilda2_ext = b12*exp(-0.5*num2_ext.^2./den2_ext.^2);

FMltilda     = FMtilda1_ext+FMtilda2_ext+FMtilda3_ext;

tan_val_incr = (0.5* tanh(1000*(5.7e-3-stretch))+0.5).*FMltilda*a_ext*kSRS.*stretch; 
tan_val_plat = (0.5* tanh(1000*(stretch-5.7e-3))+0.5).*FMltilda*a_ext*kSRS*5.7e-3; 


Fsrs = tan_val_incr + tan_val_plat; 

end

