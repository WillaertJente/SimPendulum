%% Write Results to excel 
% Write results of muscle model simulations to excel
% Jente Willaert 06-04-2021
% Gemaakt voor NCM poster

%% Import data 
map = 'C:\Users\u0125183\Box\PhD 1\Simulations Pendulum Test\Results/';
pp  = {'TD5','TD12','CP2','CP4','CP6','CP7','CP8','CP9','CP10','CP11','CP12','CP14','CP15','CP16'}; 

for j = 1:length(pp)
    for i = 1:20; 
        
        FullFileName = [map,'Result_',char(pp{j}),'_T',num2str(i),'_Fmv.mat']
        if exist(FullFileName)
            load(FullFileName)
            
            a_ext(i,j)    = sol_aext0; 
            a_flex(i,j)   = sol_a_flex; 
            kFpe_ext(i,j) = sol_kFpe_ext; 
            kFpe_flex(i,j)= sol_kFpe_flex;
            Rk(i,j)       = sol_Rk; 
            
            figure()
            plot(q_exp)
            
        end
    end
end

%% Write to excel

xlswrite('C:\Users\u0125183\Box\PhD 1\Presentaties\NCM 2021/Results_MuscleModel_Fmv.xlsx',a_ext,'a_ext')
xlswrite('C:\Users\u0125183\Box\PhD 1\Presentaties\NCM 2021/Results_MuscleModel_Fmv.xlsx',a_flex,'a_flex')
xlswrite('C:\Users\u0125183\Box\PhD 1\Presentaties\NCM 2021/Results_MuscleModel_Fmv.xlsx',kFpe_ext,'kFpe_Ext')
xlswrite('C:\Users\u0125183\Box\PhD 1\Presentaties\NCM 2021/Results_MuscleModel_Fmv.xlsx',kFpe_flex,'kFpe_Flex')
xlswrite('C:\Users\u0125183\Box\PhD 1\Presentaties\NCM 2021/Results_MuscleModel_Fmv.xlsx',Rk,'Rk')

