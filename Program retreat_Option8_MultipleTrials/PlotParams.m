function g = PlotParams(R,info)
% Plot params of tracking simulations

% Info on figure
g     = figure('Name','Optimized Parameters');
color = [0.6 0.6 0.6]; 
color_p1 = [0 95 115]./255;
color_p2 = [10 147 150]./255;

subplot(151)
bar(1,R.a_ext(1),'FaceColor',color_p1,'EdgeColor',color_p1); hold on
bar(2,R.a_ext(2),'FaceColor',color_p2,'EdgeColor',color_p2); hold on
line([0 3],[R.bounds.Ub.a R.bounds.Ub.a],'Color','k','LineWidth',1.5);hold on %Ub
line([0 3],[R.bounds.Lb.a R.bounds.Lb.a],'Color','k','LineWidth',1.5);hold on %Lb
title('a Ext'); xticks([1 2]); xticklabels({'P1','P2'})
box off;

subplot(152)
bar(1,R.a_flex(1),'FaceColor',color_p1,'EdgeColor',color_p1); hold on
bar(2,R.a_flex(2),'FaceColor',color_p2,'EdgeColor',color_p2); hold on
line([0 3],[R.bounds.Ub.a R.bounds.Ub.a],'Color','k','LineWidth',1.5);hold on %Ub
line([0 3],[R.bounds.Lb.a R.bounds.Lb.a],'Color','k','LineWidth',1.5);hold on %Lb
title('a Flex'); xticks([1 2]); xticklabels({'P1','P2'})
box off; 

subplot(153)
bar(1,R.kFpe,'FaceColor',color,'EdgeColor',color); hold on
line([0 2],[R.bounds.Ub.kFpe R.bounds.Ub.kFpe],'Color','k','LineWidth',1.5); hold on; %Ub
line([0 2],[R.bounds.Lb.kFpe R.bounds.Lb.kFpe],'Color','k','LineWidth',1.5); hold on; %Lb
line([0 2],[0.1 0.1],'Color','r','LineWidth',1.5); hold on; % nominal
ylim([0 0.2]);
title('kFpe')
box off; 

subplot(154)
bar(1,R.B,'FaceColor',color,'EdgeColor',color); hold on
line([0 2],[R.bounds.Ub.B R.bounds.Ub.B],'Color','k','LineWidth',1.5); hold on; %Ub
line([0 2],[R.bounds.Lb.B R.bounds.Lb.B],'Color','k','LineWidth',1.5); hold on; %Lb
title('B')
box off; 

subplot(155)
bar(1,R.kR(1),'FaceColor',color_p1,'EdgeColor',color_p1); hold on
bar(2,R.kR(2),'FaceColor',color_p2,'EdgeColor',color_p2); hold on
line([0 2],[R.bounds.Ub.kR R.bounds.Ub.kR],'Color','k','LineWidth',1.5); hold on % Ub
line([0 2],[R.bounds.Lb.kR R.bounds.Lb.kR],'Color','k','LineWidth',1.5); hold on % Ub
title('kR')
box off 


suptitle([info.subj ' Trial ', num2str(info.trial)])
end
