function [ xdot ] = pendulumStateDerivativeMuscleSRS_Yfb(t, x , inputdata, params )
% State derivative of muscle-driven pendulum model
% Feedback from the derivative of muscle force (Y)
% x = [q qdot lMtilde_ext Fsrs Fsrs_delayed lMtilde_flex]'
%
% Friedl De Groote - % August 24, 2017
% Jente Willaert (adaptations) - 28 may 2021 
% Friedl De Groote - further adaptations - June 7-8
global hist 

% Subject parameters 
I = params.I_OS;
m = params.mass_OS;
l = params.lc_OS;
g = 9.81; % gravity

% Inputdata = params we want to sweep 
ab_ext   = inputdata.ab_ext;         % Baseline tone extensor
ab_flex  = inputdata.ab_flex;        % Baseline tone flexor
Rk       = inputdata.Rk;             % reflex gains
kY       = inputdata.kY;
B        = inputdata.B;              % damping joint
act      = inputdata.act;
% act = interp1(inputdata.tinput, inputdata.act, t);

% moment arms and muscle tendon lengths 
m_offset = 0; 

[lMT_ext, ma_ext] = computeLengthMomentArm(params.coeff_LMT_ma_ext,x(1)+m_offset);
[lMT_flex, ma_flex] = computeLengthMomentArm(params.coeff_LMT_ma_flex,x(1)+m_offset);

% Activation 
Yfb = max(0, kY * x(7));
a_ext  = min(ab_ext + Rk * x(5) + Yfb, 1);
a_flex = ab_flex;

% state: fiber length - extensors 
[dlMtilde_ext,FT_ext,dFsrs, ~] = FiberLengthSRSOde_ext(a_ext,ab_ext,x(3)',lMT_ext, params.MTparams_ext, params.Fvparam, params.Fpparam, params.Faparam, params, x(2),x(4), inputdata); 

% state: fiber length - flexors - Note the difference - no SRS.
[dlMtilde_flex,FT_flex] = FiberLengthOde(a_flex, x(6)',lMT_flex, params.MTparams_flex, params.Fvparam, params.Fpparam, params.Faparam);

% xdot = derivative of x (= q qdot lMtilda_ext Fsrs Fsrs_delayed lMtilda_flex)
xdot    = zeros(6,1);
xdot(1) = x(2);                                                                     % xdot
xdot(2) = 1/I * ((-m*g*l*cos(x(1)))- FT_ext.*ma_ext - FT_flex.*ma_flex - B*x(2));    % xdd
xdot(3) = dlMtilde_ext';                                                            % Derivative of fiber length (extensor)-  From FiberLengthSRSOde 
xdot(4) = dFsrs;                                                                    %  Derivative of SRS force - From FiberLengthSRSode
xdot(5) = (x(4) - x(5))/inputdata.tau_d;                                            % Derivative of delayed SRS force.
xdot(6) = dlMtilde_flex'; % Derivative of fiber length (Flexor) -  From FiberLengthOde.
xdot(7) = (dFsrs - x(7))/inputdata.tau_d;
end

