function [ xdot_ext] = isometricDerivative_polynomial_ext(t, x , params,inputdata )
% State derivative of muscle-driven pendulum model
m_offset = 0; 
theta    = (-120:0.01:10)' *pi/180;

lMT_ext = params.coeff_LMT_ma_ext(1) + params.coeff_LMT_ma_ext(2)*(x+m_offset) + params.coeff_LMT_ma_ext(3)*(x+m_offset).^2 + params.coeff_LMT_ma_ext(4)*(x+m_offset).^3 + params.coeff_LMT_ma_ext(5)*(x+m_offset).^4 + params.coeff_LMT_ma_ext(6)*(x+m_offset).^5;
ma_ext  = -params.coeff_LMT_ma_ext(2) + -params.coeff_LMT_ma_ext(3)*(x+m_offset) + -params.coeff_LMT_ma_ext(4)*(x+m_offset).^2 + -params.coeff_LMT_ma_ext(5)*(x+m_offset).^3 + -params.coeff_LMT_ma_ext(6)*(x+m_offset).^4;

lMT_flex = params.coeff_LMT_ma_flex(1) + params.coeff_LMT_ma_flex(2)*(x+m_offset) + params.coeff_LMT_ma_flex(3)*(x+m_offset).^2 + params.coeff_LMT_ma_flex(4)*(x+m_offset).^3 + params.coeff_LMT_ma_flex(5)*(x+m_offset).^4 + params.coeff_LMT_ma_flex(6)*(x+m_offset).^5;
ma_flex  = -params.coeff_LMT_ma_flex(2) + -params.coeff_LMT_ma_flex(3)*(x+m_offset) + -params.coeff_LMT_ma_flex(4)*(x+m_offset).^2 + -params.coeff_LMT_ma_flex(5)*(x+m_offset).^3 + -params.coeff_LMT_ma_flex(6)*(x+m_offset).^4;

% LMT_ext_start = interp1(theta, lMT_ext, inputdata.q0);
% ma_ext_start  = interp1(theta, ma_ext, inputdata.q0);
% LMT_flex_start= interp1(theta, lMT_flex, inputdata.q0);
% ma_flex_start = interp1(theta, ma_flex, inputdata.q0); 

vMT_ext = ma_ext * 0;               % ?
vMT_flex = ma_flex *0;

a.ext  = inputdata.ab_ext;
a.flex = inputdata.ab_flex; 

% state: fiber length
[dx_ext,FT_ext, dx_flex, FT_flex] = FiberLengthOde_2muscles(a,x',lMT_ext, lMT_flex, params.MTparams_ext, params.MTparams_flex , params.Fvparam, params.Fpparam, params.Faparam);

% state: fiber force
% [dx, FT] = TendonForceOde(a, x', LMT, vMT, params.MTparams, params.Fvparam, params.Fpparam, params.Faparam);

xdot_ext  = dx_ext';

end

