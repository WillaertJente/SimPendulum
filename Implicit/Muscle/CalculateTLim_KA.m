function [Tlim] = CalculateTLim_KA(x, Kr1, Kr2, klim, offset)
%Function to calculate T limit
%   Tlim = passive force to prevent joint from going outside their range of
%   motion
%   TLim ifv knee hoek ipv hoek tov bovenbeen

offset   = offset*pi/180; 
m_offset = mean(offset);

knee_r_range = [Kr1 Kr2]*pi/180;
theta_ref    = (knee_r_range(2) - knee_r_range(1))/2;
q            = (x+m_offset) - (knee_r_range(1) + knee_r_range(2))/2;

Tlim         = -exp(klim * (q - theta_ref)) + exp(klim * (-q - theta_ref));


end

