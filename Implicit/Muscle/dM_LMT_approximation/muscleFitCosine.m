function  x  = muscleFitCosine( angle, dM, LMT )

% dM = a + b*cos(c*theta + d) 
% lMT = - (a*theta + b/c*sin(c*theta + d) + h

options = optimoptions('fminunc','Algorithm','quasi-newton');
options.MaxFunctionEvaluations = 10000;
f = @(x)getFitCost(x,pi/180*angle, dM, LMT )
[x,fval] = fminunc(f,-ones(5,1),options);




end


