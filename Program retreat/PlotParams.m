function g = PlotParams(R,info)
% Plot params of tracking simulations

% Info on figure
g     = figure('Name','Optimized Parameters');
color = [0.6 0.6 0.6]; 

subplot(121)
bar(1,R.a,'FaceColor',color,'EdgeColor',color); hold on
line([0 2],[R.bounds.Ub.a R.bounds.Ub.a],'Color','k','LineWidth',1.5);hold on %Ub
line([0 2],[R.bounds.Lb.a R.bounds.Lb.a],'Color','k','LineWidth',1.5);hold on %Lb
title('a')
%title([info.subj ' Trial ', num2str(info.trial)])
box off; 

subplot(122)
bar(1,R.kFpe,'FaceColor',color,'EdgeColor',color); hold on
line([0 2],[R.bounds.Ub.kFpe R.bounds.Ub.kFpe],'Color','k','LineWidth',1.5); hold on; %Ub
line([0 2],[R.bounds.Lb.kFpe R.bounds.Lb.kFpe],'Color','k','LineWidth',1.5); hold on; %Lb
line([0 2],[0.1 0.1],'Color','r','LineWidth',1.5); hold on; % nominal
title('kFpe')
box off; 

suptitle([info.subj ' Trial ', num2str(info.trial)])
end
