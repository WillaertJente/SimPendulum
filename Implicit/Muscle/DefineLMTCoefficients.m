function [coeff_LMT_ma] = DefineLMTCoefficients(map_MA, sn)
%  Calculate LMT and ma coefficients

% Theta
%theta    = (-120:0.01:10)' *pi/180; theta_sq = theta.^2; theta_th = theta.^3; theta_fo = theta.^4;
theta    = (-140:0.01:20)' *pi/180; theta_sq = theta.^2; theta_th = theta.^3; theta_fo = theta.^4; % for cp4 and cp 8 
% Name MA files 
name_LMT = '_MuscleAnalysis_Length.sto';
name_ma  = '_MuscleAnalysis_MomentArm_knee_angle_l.sto';        % CP4 TD 5
%name_ma  = '_MuscleAnalysis_MomentArm_knee_angle_r.sto'   ;     % CP 8

% Data LMT
% LMT_dat  = importdata([map_MA,char(sn),'_',name_LMT]);
LMT_dat  = importdata([map_MA,name_LMT]);   
col_RF   = find(strcmp('rect_fem_l',LMT_dat.colheaders));       % CP 4 TD 5
%col_RF   = find(strcmp('rect_fem_r',LMT_dat.colheaders));      % CP 8
LMT      = LMT_dat.data(:,col_RF);

% Data ma
% ma_dat   = importdata([map_MA,char(sn),'_',name_ma]);
ma_dat   = importdata([map_MA,name_ma]);
col_RF   = find(strcmp('rect_fem_l',ma_dat.colheaders));        % CP 4 TD 5
% col_RF   = find(strcmp('rect_fem_r',ma_dat.colheaders));      % CP 8
ma       = ma_dat.data(:,col_RF);

% Define LMT coefficients together with ma coefficients
LMT_ma_dat = [LMT; -ma];
one      = ones(length(LMT),1);
zero     = zeros(length(ma),1); 

A        = [one  theta theta_sq   theta_th ];
B        = [zero one   2*theta    3*theta_sq ];
C        = [A; B];

coeff_LMT_ma = C\LMT_ma_dat;

res_LMT = sqrt(mean((A*coeff_LMT_ma - LMT).^2))
res_ma  = sqrt(mean((-B*coeff_LMT_ma - ma).^2))

res1     = A*coeff_LMT_ma - LMT;
res2     = -B*coeff_LMT_ma - ma;
plot(theta,res1)
hold on
plot(theta,res2)

% lMT = coeff_LMT_ma(1)  + coeff_LMT_ma(2)*x + coeff_LMT_ma(3)*x.^2 + coeff_LMT_ma(4)*x.^3;
% ma  = -coeff_LMT_ma(2) + -coeff_LMT_ma(3)*x + -coeff_LMT_ma(4)*x.^2;
end