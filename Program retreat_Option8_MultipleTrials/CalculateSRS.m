function [Fsrs, dFsrs2dt] = CalculateSRS(info, data_exp, lMtilda, a_ext, Fsrs2, FMltilda)
% Calculate Fsrs and Derivative of FSRS 

% Input 
kSRS    = info.kSRS; 
N_1     = data_exp.N_1; 
N       = data_exp.Nspline;

% Stretch 
stretch = lMtilda(1,1:N_1)-lMtilda(1,1); 

% Fsrs 1
Fsrs1   = stretch.*FMltilda(1,1:N_1)*a_ext*kSRS; 

% FSRS
Fsrs    = [Fsrs1 Fsrs2]; 

%dFsrsdt
dFsrs2dt = Fsrs/0.05;  

end

