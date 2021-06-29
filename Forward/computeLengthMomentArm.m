function [lMT, ma] = computeLengthMomentArm(c,x)
% Friedl De Groote
% June 7, 2021
%
% Inputs
% c     vector of length n with coefficients
% x     joint angle
%
% Outputs
% Muscle tendon length and moment arm
% lMT = c1 * x^0 + c2 * x^1 + c3 * x^2 + ... + cn * x^(n-1)
% ma = c2 * x^0 + 2 * c3 * x^1 + 3 * c4 * x^2 + ... + (n-1) * cn * x^(n-2)

lMT = c(1)*ones(size(x));
ma = 0*ones(size(x));

for i = 1:length(c)-1
   lMT = lMT + c(i+1)*x.^i;
   ma = ma + i*c(i+1)*x.^(i-1);
end