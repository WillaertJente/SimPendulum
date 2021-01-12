function [ xdot] = isometricDerivative(t, x , params )
% State derivative of muscle-driven pendulum model
% x = [q qdot dF_RF dF_VAS]'
%
% Friedl De Groote
% August 24, 2017

LMT = interp1(params.LMT_qknee, params.LMT, params.q0);
ma = interp1(params.ma_qknee, params.ma, params.q0);
vMT = ma * 0;

a = params.ab;

% state: fiber length
[dx,FT] = FiberLengthOde(a,x',LMT, params.MTparams, params.Fvparam, params.Fpparam, params.Faparam);

% state: fiber force
% [dx, FT] = TendonForceOde(a, x', LMT, vMT, params.MTparams, params.Fvparam, params.Fpparam, params.Faparam);

xdot = dx';

end

