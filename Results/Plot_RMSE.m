%% Plot RMSE 

RMSE.TD5 = [1.065701 16.90798]; 
RMSE.CP2 = [3.437747 7.081758 2.463719 3.220023]; 
RMSE.CP8 = [1.827735 1.105809 1.644389 1.9022];
RMSE.CP9 = [4.222699 2.887707 6.10773];
RMSE.CP10= [4.291454 2.939273 3.517961];
RMSE.CP11= [3.701307 2.108485 3.300237];
RMSE.SDR = [2.801764 1.799087 5.723848 2.320479 3.323155];

mean = 3.164969; 

figure()
subplot(221)
bar(1,mean,'FaceColor',[0.7 0.7 0.7],'EdgeColor',[0.7 0.7 0.7]); hold on
scatter(ones(1,2)*0.9, RMSE.TD5,'MarkerEdgeColor',[0 119 182]./255,'MarkerFaceColor',[0 119 182]./255) % TD5
scatter(ones(1,4)*0.95, RMSE.CP2,'MarkerEdgeColor',[255 183 3]./255,'MarkerFaceColor',[255 183 3]./255) % CP2
scatter(ones(1,4)*1, RMSE.CP8,'MarkerEdgeColor',[251 133 0]./255,'MarkerFaceColor',[251 133 0]./255) % CP8
scatter(ones(1,3)*1.05, RMSE.CP9,'MarkerEdgeColor',[229 107 111]./255,'MarkerFaceColor',[229 107 111]./255) % CP9
scatter(ones(1,3)*1.1, RMSE.CP10,'MarkerEdgeColor',[214 40 40]./255,'MarkerFaceColor',[214 40 40]./255) % CP10
scatter(ones(1,3)*1.15, RMSE.CP11,'MarkerEdgeColor',[252 191 73]./255,'MarkerFaceColor',[252 191 73]./255) % CP11
scatter(ones(1,5)*0.85, RMSE.SDR,'MarkerEdgeColor',[173 193 120]./255,'MarkerFaceColor',[173 193 120]./255) % CP11
box off; title('RMSE');
