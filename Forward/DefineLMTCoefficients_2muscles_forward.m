function [coeff_LMT_ma_ext, coeff_LMT_ma_flex] = DefineLMTCoefficients_2muscles_forward(map_MA, sn, vis)
%  Calculate LMT and ma coefficients

% Theta of fake motion file 
if strcmp('CP6',sn) % TD12
    theta = (-140:0.5:20)' *pi/180; theta_sq = theta.^2; theta_th = theta.^3; theta_fo = theta.^4; theta_fi = theta.^5;
elseif strcmp('TD5',sn) 
    theta = (-120:0.01:10)' *pi/180; theta_sq = theta.^2; theta_th = theta.^3; theta_fo = theta.^4; theta_fi = theta.^5;
else
    theta = (-140:0.01:20)' *pi/180; theta_sq = theta.^2; theta_th = theta.^3; theta_fo = theta.^4; theta_fi = theta.^5;    
end

% Name MA files
name_LMT = '_MuscleAnalysis_Length.sto';
if strcmp('CP6',sn)
    name_ma  = '_MuscleAnalysis_MomentArm_knee_angle_r.sto'; 
else
    name_ma = '_MuscleAnalysis_MomentArm_knee_angle_l.sto'; 
end

% Data LMT
LMT_dat  = importdata([map_MA,name_LMT]);   
if strcmp('CP6',sn)
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
if strcmp('CP6',sn)
    col_RF  = find(strcmp('rect_fem_r',ma_dat.colheaders));
    col_BF  = find(strcmp('bifemlh_r',ma_dat.colheaders));
else
    col_RF  = find(strcmp('rect_fem_l',ma_dat.colheaders));
    col_BF  = find(strcmp('bifemlh_l',ma_dat.colheaders));
end
ma(:,1)     = ma_dat.data(:,col_RF);
ma(:,2)     = ma_dat.data(:,col_BF);

% Define LMT coefficients together with ma coefficients
one      = ones(length(LMT),1);
zero     = zeros(length(ma),1); 

A        = [one  theta theta_sq   theta_th    theta_fo    theta_fi];
B        = [zero one   2*theta    3*theta_sq  4*theta_th  5*theta_fo ];
C        = [A; B];

coeff_LMT_ma = zeros(6,2);

for i = 1:2
    LMT_ma_dat = [LMT(:,i); -ma(:,i)];
    coeff_LMT_ma(:,i) = C\LMT_ma_dat;
    
    res_LMT = sqrt(mean((A*coeff_LMT_ma(:,i) - LMT(:,i)).^2))
    res_ma  = sqrt(mean((-B*coeff_LMT_ma(:,i) - ma(:,i)).^2))
    
    if vis == 1
        figure()
        subplot(211)
        plot(theta,LMT(:,i),'Color',[0.5 0.5 0.5],'LineWidth', 3); hold on;
        plot(theta,A*coeff_LMT_ma(:,i),'k','LineWidth', 1);
        ylabel('lMT')
        title('extensor')
        subplot(212)
        plot(theta,-ma(:,i),'Color',[0.5 0.5 0.5],'LineWidth', 3); hold on;
        plot(theta,B*coeff_LMT_ma(:,i),'k','LineWidth', 1);
        xlabel('theta')
        ylabel('moment arm')
        legend('input','fit')
    end
    
end
coeff_LMT_ma_ext = coeff_LMT_ma(:,1);
coeff_LMT_ma_flex = coeff_LMT_ma(:,2);
