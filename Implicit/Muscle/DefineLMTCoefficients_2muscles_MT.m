function [coeff_LMT_ma_ext, coeff_LMT_ma_flex] = DefineLMTCoefficients_2muscles_MT(map_MA, sn)
%  Calculate LMT and ma coefficients

% Theta
if strcmp('CP8',sn) % TD12
    theta = (-140:0.5:20)' *pi/180; theta_sq = theta.^2; theta_th = theta.^3; theta_fo = theta.^4; theta_fi = theta.^5;
elseif strcmp('TD5',sn) 
    theta = (-120:0.01:10)' *pi/180; theta_sq = theta.^2; theta_th = theta.^3; theta_fo = theta.^4; theta_fi = theta.^5;
else
    theta = (-140:0.01:20)' *pi/180; theta_sq = theta.^2; theta_th = theta.^3; theta_fo = theta.^4; theta_fi = theta.^5;    
end

% Name MA files
name_LMT = '_MuscleAnalysis_Length.sto';
if strcmp('CP8',sn)
    name_ma  = '_MuscleAnalysis_MomentArm_knee_angle_r.sto'; 
else
    name_ma = '_MuscleAnalysis_MomentArm_knee_angle_l.sto'; 
end

% Data LMT
% LMT_dat  = importdata([map_MA,char(sn),'_',name_LMT]);
LMT_dat  = importdata([map_MA,name_LMT]);   
if strcmp('CP8',sn)
    col_RF   = find(strcmp('rect_fem_r',LMT_dat.colheaders));
    col_BF   = find(strcmp('bifemlh_r',LMT_dat.colheaders));
else
    col_RF   = find(strcmp('rect_fem_l',LMT_dat.colheaders));
    col_BF   = find(strcmp('bifemlh_l',LMT_dat.colheaders));
end
LMT(:,1)     = LMT_dat.data(:,col_RF);
LMT(:,2)     = LMT_dat.data(:,col_BF);

% Data ma
ma_dat   = importdata([map_MA,name_ma]);
if strcmp('CP8',sn)
    col_RF  = find(strcmp('rect_fem_r',ma_dat.colheaders));
    col_BF  = find(strcmp('bifemlh_r',ma_dat.colheaders));
else
    col_RF  = find(strcmp('rect_fem_l',ma_dat.colheaders));
    col_BF  = find(strcmp('bifemlh_l',ma_dat.colheaders));
end
ma(:,1)     = ma_dat.data(:,col_RF);
ma(:,2)     = ma_dat.data(:,col_BF);

% Define LMT coefficients together with ma coefficients    for extensors 
LMT_ma_dat_ext = [LMT(:,1); -ma(:,1)];
one      = ones(length(LMT),1);
zero     = zeros(length(ma),1); 

A        = [one  theta theta_sq   theta_th    theta_fo    theta_fi];
B        = [zero one   2*theta    3*theta_sq  4*theta_th  5*theta_fo ];
C        = [A; B];

coeff_LMT_ma_ext = C\LMT_ma_dat_ext;

res_LMT_ext = sqrt(mean((A*coeff_LMT_ma_ext - LMT(:,1)).^2))
res_ma_ext  = sqrt(mean((-B*coeff_LMT_ma_ext - ma(:,1)).^2))

res1_ext     = A*coeff_LMT_ma_ext - LMT(:,1);
res2_ext     = -B*coeff_LMT_ma_ext - ma(:,1);

% Define LMT coefficients together with ma coefficients    for flexors 
LMT_ma_dat_flex = [LMT(:,2); -ma(:,2)];
one      = ones(length(LMT),1);
zero     = zeros(length(ma),1); 

A        = [one  theta theta_sq   theta_th    theta_fo   theta_fi];
B        = [zero one   2*theta    3*theta_sq  4*theta_th 5*theta_fo];
C        = [A; B];

coeff_LMT_ma_flex = C\LMT_ma_dat_flex;

res_LMT_flex = sqrt(mean((A*coeff_LMT_ma_flex - LMT(:,2)).^2))
res_ma_flex  = sqrt(mean((-B*coeff_LMT_ma_flex - ma(:,2)).^2))

res1_flex     = A*coeff_LMT_ma_flex - LMT(:,2);
res2_flex     = -B*coeff_LMT_ma_flex - ma(:,2);

% plot(theta,res1)
% hold on
% plot(theta,res2)

% lMT = coeff_LMT_ma(1)  + coeff_LMT_ma(2)*x + coeff_LMT_ma(3)*x.^2 + coeff_LMT_ma(4)*x.^3;
% ma  = -coeff_LMT_ma(2) + -coeff_LMT_ma(3)*x + -coeff_LMT_ma(4)*x.^2;
end