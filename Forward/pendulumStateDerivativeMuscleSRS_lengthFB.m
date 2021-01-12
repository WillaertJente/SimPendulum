function [ xdot ] = pendulumStateDerivativeMuscle(t, x, Z, dZ)
% State derivative of muscle-driven pendulum model
% x = [q qdot dF_RF dF_VAS]'
%
% Friedl De Groote
% August 24, 2017

global params

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

ma = interp1(params.ma_qknee, params.ma, x(1));
LMT = interp1(params.LMT_qknee, params.LMT, x(1));
vMT = ma * x(2);

e = params.ab + max(0, params.kl * (Z(5)-params.lMiso) + params.kv * dZ(5));
dadt = ActivationOde(e,x(4));

% state: fiber length
[dx,FT,Fsrs] = FiberLengthSRSOde(x(4)', x(5)', LMT, params.MTparams, params.Fvparam, params.Fpparam, params.Faparam, params, x(2));

% state: fiber force
% [dx, FT] = TendonForceOde(a, x', LMT, vMT, params.MTparams, params.Fvparam, params.Fpparam, params.Faparam);

% SRS
xdot = zeros(3,1);
xdot(1) = x(2);
xdot(2) = 1/I *(-m*g*l*cos(x(1)) + Tlimit + x(3) + FT(1)*ma(1)); % + FT(2)*ma(2));
xdot(3) = params.kM*(-x(2) - x(3)/(params.d));
xdot(4) = dadt;
xdot(5) = dx';

end

