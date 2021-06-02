%% Plot difference between multiple trials and single trials 
% Jente Willaert - 02062021

% Load data 
map_result = 'C:\Users\u0125183\Box\PhD 1\Simulations Pendulum Test\Results/Result_CP7_'; 

name.s = 'T5_ScaledTorso.mat';  % data single trial
name.m = 'MT8.mat';             % data multiple trials
nm     = 2                   % trialnummer in multiple oplossing
name.m2= 'MT10.mat'; 
nm2    = 1

data.s = load([map_result, name.s]);
data.m = load([map_result, name.m]);
data.m2 = load([map_result, name.m2]);

% Plot results 
figure()
plot(data.s.tvect, data.s.q_exp,'k','LineWidth',1.5); hold on
plot(data.s.tvect, data.s.sol_x,'Color',[69,123,157]./255,'LineWidth',1.5); hold on
plot(data.s.tvect, data.m.sol_x(1:length(data.s.tvect)),'Color',[214,40,40]./255,'LineWidth',1.5); hold on
plot(data.s.tvect, data.m2.sol_x(1:length(data.s.tvect)),'LineWidth',1.5); hold on; 
legend({'Exp','Single','Multiple','Multiple2'})

figure()
subplot(251)
bar(1,data.s.sol_aext0,'EdgeColor',[69,123,157]./255,'FaceColor',[69,123,157]./255); hold on
bar(2,data.m.sol_aext0(nm),'EdgeColor',[214,40,40]./255,'FaceColor',[214,40,40]./255); hold on
bar(3,data.m2.sol_aext0(nm2)); hold on
box off; legend({'s','m','m2'}); title('a ext'); 

subplot(252)
bar(1,data.s.sol_a_flex,'EdgeColor',[69,123,157]./255,'FaceColor',[69,123,157]./255); hold on
bar(2,data.m.sol_a_flex(nm),'EdgeColor',[214,40,40]./255,'FaceColor',[214,40,40]./255); hold on
bar(3,data.m2.sol_a_flex(nm2)); hold on
box off;  title('a flex'); 

subplot(253)
bar(1,data.s.sol_Rk,'EdgeColor',[69,123,157]./255,'FaceColor',[69,123,157]./255); hold on
bar(2,data.m.sol_Rk(nm),'EdgeColor',[214,40,40]./255,'FaceColor',[214,40,40]./255); hold on
bar(3,data.m2.sol_Rk(nm2)); hold on
box off;  title('Rk'); 

subplot(254)
bar(1,data.s.sol_kFpe_ext,'EdgeColor',[69,123,157]./255,'FaceColor',[69,123,157]./255); hold on
bar(2,data.m.sol_kFpe_ext,'EdgeColor',[214,40,40]./255,'FaceColor',[214,40,40]./255); hold on
bar(3,data.m2.sol_kFpe_ext(nm2)); hold on
box off; title('kFpe ext'); 

subplot(255)
bar(1,data.s.sol_kFpe_flex,'EdgeColor',[69,123,157]./255,'FaceColor',[69,123,157]./255); hold on
bar(2,data.m.sol_kFpe_flex,'EdgeColor',[214,40,40]./255,'FaceColor',[214,40,40]./255); hold on
bar(3,data.m2.sol_kFpe_flex(nm2)); hold on
box off; title('kFpe flex'); 
