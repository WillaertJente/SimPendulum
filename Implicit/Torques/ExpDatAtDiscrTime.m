function [q_exp, qdot_exp, N, tvect, dt] = ExpDatAtDiscrTime(t_span,t_exp,q_exp_r)  
%Interpolate data at discretized time using spline
%   1. Get time vector for simulation (tvect)
%   2. Interpolate experimental data 

dt    = 0.005;                              % Time step of simulation
tvect = t_span(1):dt:t_span(end);           % Time vector for simulation
N     = length(tvect) ;                      % Number of steps in simulation

qspline = spline(t_exp,q_exp_r);            % Spline fit
[q_exp,qdot_exp] = SplineEval_ppuval(qspline,tvect,1); % Get angles and velocities
end

