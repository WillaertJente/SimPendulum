function dx = dldtIsometric(t, x , params,inputdata, mus )
% This function returns the derivative of muscle length (x) at a constant knee
% angle (inputdata.q0).
m_offset = 0; 

eval(['coeff = params.coeff_LMT_ma_', mus,';']) 
[lMT, ~] = computeLengthMomentArm(coeff,inputdata.q0+m_offset);

eval(['a  = inputdata.ab_', mus,';'])
eval(['MTparams  = params.MTparams_', mus,';'])

% state: fiber length
[dx,~] = FiberLengthOde(a,x,lMT, MTparams, params.Fvparam, params.Fpparam, params.Faparam);

end

