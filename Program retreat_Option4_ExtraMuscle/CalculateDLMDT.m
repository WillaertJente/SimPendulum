function [dlMdt] = CalculateDLMDT(vMtilda, params_OS)
%Calculate derivative of muscle length

lMo        = params_OS.MT(2,:); 
vMtildamax = params_OS.MT(5,:);
dlMdt_ext  = vMtilda.* vMtildamax(1)./ lMo(1);  % Extensor - change vMtilda
dlMdt_flex = vMtilda.* vMtildamax(2)./ lMo(2);  % Flexor - change vMtilda
dlMdt      = dlMdt_ext; 
end

