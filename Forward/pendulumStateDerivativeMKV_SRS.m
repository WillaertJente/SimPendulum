function [ xdot ] = pendulumStateDerivative(t, x , params )
% State derivative of pendulum model
% x = [q qdot Td TSRS]'
global hist

I = params.I;
m = params.mass;
l = params.lc;
g = 9.81; % gravity

% Passive force to prevent joints from going outside their range of motion.
knee_r_range = params.knee_r_range;
klim = params.klim;
theta_ref = (knee_r_range(2) - knee_r_range(1))/2;
q = x(1) - (knee_r_range(1) + knee_r_range(2))/2;
Tlimit = -exp(klim*(q-theta_ref)) + exp(klim*(-q-theta_ref));

% Maxwell model
% Td = - params.d * x(2);
Td = 0;

% SRS
xdot = zeros(4,1);
xdot(1) = x(2);
xdot(2) = 1/I *(-m*g*l*cos(x(1)) + Tlimit + x(3) + x(4) + Td + params.Tb);
xdot(3) = params.kM*(-x(2) - x(3)/(params.d));

kSRS = params.kSRS * params.Tb; 
delta_theta = x(1) - params.theta0;
% if abs(delta_theta) < params.delta_theta_crit && x(2) <= 0 && hist == 1
%     xdot(4) = 1/0.001 * (-x(4) - kSRS * delta_theta);
% elseif x(2) <= 0 && hist == 1
%     xdot(4) = 1/0.001 * (-x(4) + kSRS * params.delta_theta_crit);
% else 
%     xdot(4) = 1/params.tauSRS * (-x(4));
%     hist = 0;
% end
% Smooth version
% if x(2) <= 0 & hist == 1
%     s_stcrit = 1/2*(erf(100*(-delta_theta + params.delta_theta_crit)) + 1);
%     s_gtcrit = 1/2*(erf(100*(delta_theta - params.delta_theta_crit)) + 1);    
%     xdot(4) = 1/0.001 * ( s_stcrit*(-x(4) - kSRS * delta_theta) + s_gtcrit*(-x(4) + kSRS * params.delta_theta_crit));
% else 
%     xdot(4) = 1/params.tauSRS * (-x(4));
%     hist = 0;
% end
s_velneg = 1/2*(erf(-10*x(2)) + 1);
if hist == 1 & s_velneg > 0.001
    s_stcrit = 1/2*(erf(50*(delta_theta + params.delta_theta_crit)) + 1);
    s_gtcrit = 1/2*(erf(50*(-delta_theta- params.delta_theta_crit)) + 1);    
    xdot_SRS = 1/0.001 * ( s_stcrit*(-x(4) - kSRS * delta_theta) + s_gtcrit*(-x(4) + kSRS * params.delta_theta_crit));
    s_velpos = 1/2*(erf(10*x(2)) + 1);
    xdot(4) = s_velneg * xdot_SRS + s_velpos * 1/params.tauSRS * (-x(4));   

else
    xdot(4) = 1/params.tauSRS * (-x(4));
    hist = 0;
end

end

