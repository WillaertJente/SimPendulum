function [ xdot ] = pendulumStateDerivativeMuscleSRS_howtoproceed(t, x , inputdata, params )
% State derivative of muscle-driven pendulum model
% x = [q qdot lMtilde_ext Fsrs Fsrs_delayed lMtilde_flex]'
%
% Friedl De Groote - % August 24, 2017
% Jente Willaert (adaptations) - 28 may 2021 
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
B        = inputdata.B;              % damping joint
kFpe_ext = inputdata.kFpe_ext;       % Passive force extensor
kFpe_flex= inputdata.kFpe_flex;      % Passive force flexor 
act      = inputdata.act;

% moment arms and muscle tendon lengths 
m_offset = 0; 
theta    = (-120:0.01:10)' *pi/180;

lMT_ext = params.coeff_LMT_ma_ext(1) + params.coeff_LMT_ma_ext(2)*(x(1)+m_offset) + params.coeff_LMT_ma_ext(3)*(x(1)+m_offset).^2 + params.coeff_LMT_ma_ext(4)*(x(1)+m_offset).^3 + params.coeff_LMT_ma_ext(5)*(x(1)+m_offset).^4 + params.coeff_LMT_ma_ext(6)*(x(1)+m_offset).^5;
ma_ext  = -params.coeff_LMT_ma_ext(2) + -params.coeff_LMT_ma_ext(3)*(x(1)+m_offset) + -params.coeff_LMT_ma_ext(4)*(x(1)+m_offset).^2 + -params.coeff_LMT_ma_ext(5)*(x(1)+m_offset).^3 + -params.coeff_LMT_ma_ext(6)*(x(1)+m_offset).^4;

lMT_flex = params.coeff_LMT_ma_flex(1) + params.coeff_LMT_ma_flex(2)*(x(1)+m_offset) + params.coeff_LMT_ma_flex(3)*(x(1)+m_offset).^2 + params.coeff_LMT_ma_flex(4)*(x(1)+m_offset).^3 + params.coeff_LMT_ma_flex(5)*(x(1)+m_offset).^4 + params.coeff_LMT_ma_flex(6)*(x(1)+m_offset).^5;
ma_flex  = -params.coeff_LMT_ma_flex(2) + -params.coeff_LMT_ma_flex(3)*(x(1)+m_offset) + -params.coeff_LMT_ma_flex(4)*(x(1)+m_offset).^2 + -params.coeff_LMT_ma_flex(5)*(x(1)+m_offset).^3 + -params.coeff_LMT_ma_flex(6)*(x(1)+m_offset).^4;

vMT_ext  = ma_ext * x(2);
vMT_flex = ma_flex * x(2);

% Activation 
a_ext  = ab_ext + Rk * x(5);
a_flex = ab_flex;

% state: fiber length - extensors 
[dlMtilde_ext,FT_ext,dFsrs, Fpe] = FiberLengthSRSOde_ext(a_ext,x(3)',lMT_ext, params.MTparams_ext, params.Fvparam, params.Fpparam, params.Faparam, params, x(2),x(4),kFpe_ext, inputdata); % kFpe toeveogen

% state: fiber length - flexors - Note the difference - no SRS.
[dlMtilde_flex,FT_flex] = FiberLengthOde_flex(a_flex, x(6)',lMT_flex, params.MTparams_flex, params.Fvparam, params.Fpparam, params.Faparam, kFpe_flex);

% xdot = derivative of x (= q qdot lMtilda_ext Fsrs Fsrs_delayed lMtilda_flex)
xdot    = zeros(6,1);
xdot(1) = x(2);                                                                     % xdot
xdot(2) = 1/I * ((-m*g*l*cos(x(1)))+ FT_ext.*ma_ext + FT_flex.*ma_flex - B*x(2) );    % xdd
xdot(3) = dlMtilde_ext';                                                            % Derivative of fiber length (extensor)-  From FiberLengthSRSOde 
xdot(4) = dFsrs;                                                                    %  Derivative of SRS force - From FiberLengthSRSode
xdot(5) = (x(4) - x(5))/inputdata.tau_d;                                            % Derivative of delayed SRS force.
xdot(6) = dlMtilde_flex';                                                           % Derivative of fiber length (Flexor) -  From FiberLengthOde.
end

