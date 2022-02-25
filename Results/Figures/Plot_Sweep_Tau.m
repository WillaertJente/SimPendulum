%% Plot tau sweep 
color = [0 0 0; 216 227 233; 189 208 219; 150 180 197; 110 151 171; 80 121 145; 65 99 118; 51 77 92;]./255
for i = 2:8 
    Name = ['C:\Users\u0125183\Documents\MATLAB\SimPendulum - programming retreat\Results/CP8_T2_Opt7_Reflexes_Tau0.0',num2str(i),'.mat']
    load(Name); 
    
    if i == 2
    figure(1)
    plot(R.exp.qspline,'k','LineWidth',1.5); hold on
    plot(R.x,'Color',color(i,:),'LineWidth',1.5); hold on
    else
        figure(1)
        plot(R.x,'Color',color(i,:),'LineWidth',1.5); hold on
    end
    figure(1)
    legend({'Exp','Tau 0.02','0.03','0.04','0.05','0.06','0.07','0.08'}); 
    box off; title('CP8 T2')
    
    figure(2)
    subplot(2,2,1)
    bar(i-1,R.a_ext,'FaceColor',color(i,:),'EdgeColor',color(i,:)); hold on
    box off; ylim([0 0.03]); title(' a ext')
    subplot(2,2,2)
    bar(i-1,R.a_flex,'FaceColor',color(i,:),'EdgeColor',color(i,:)); hold on
    box off; ylim([0 0.03]); title('a flex')
    subplot(2,2,3)
    bar(i-1,R.kFpe,'FaceColor',color(i,:),'EdgeColor',color(i,:)); hold on
    box off; ylim([0 0.21]); title('kFpe')
    subplot(2,2,4)
    bar(i-1,R.kR,'FaceColor',color(i,:),'EdgeColor',color(i,:)); hold on
    box off; ylim([0 10]); title('kR')
    
    
end
figure(2)
suptitle('CP8 T2')
