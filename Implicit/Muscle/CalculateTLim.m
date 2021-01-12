function [Tlim] = CalculateTLim(x, Kr1, Kr2, klim)
%Function to calculate T limit
%   Tlim = passive force to prevent joint from going outside their range of
%   motion

knee_r_range = [Kr1 Kr2]*pi/180;
theta_ref    = (knee_r_range(2) - knee_r_range(1))/2;
q            = x - (knee_r_range(1) + knee_r_range(2))/2;

Tlim         = -exp(klim * (q - theta_ref)) + exp(klim * (-q - theta_ref));


end

