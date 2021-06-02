function [q_exp_r, t_exp, t_span, on_srs] = LoadExpData_forward(name, trial, bool, params, path) 
%Load experimental data 
%   1. Create name of trial you want to load
%   2. Load data in structure (1 = time, 2 = knee angle (degrees) with rest
%   on 0
%   3. Export experimental angles in radians (q_exp_r)
%   4. Export experimental time (t_exp)

% s.tot   = [char(name),'_T',num2str(trial(j))];      % Name of data
s.tot   = ['Trial',num2str(trial)];
q_t_exp = load([path,'BK_',s.tot,'.mat']);           % Experimental data (time, angle in degrees, rust op 0)

q_exp_r = (q_t_exp.data(:,2)-90)*pi/180;
t_exp   = q_t_exp.data(:,1); 
t_span  = [t_exp(1) t_exp(end)];

if bool == 1
    plot(t_exp,q_exp_r,'k','LineWidth',1.5)
end

if trial < params.Nmr
    on_srs = 1;
else
    on_srs = 0; 
end

end

