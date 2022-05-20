%% Figure to investigate the effect of different IG 

% Plot trajectories to see whether there is a difference 
%%
subj  = 'TD5';
trial = 1;
s  = load([subj,'_T',num2str(trial),'_Opt7_IG_standard.mat']);
v1 = load([subj,'_T',num2str(trial),'_Opt7_IG_aflex_v1.mat']);
v2 = load([subj,'_T',num2str(trial),'_Opt7_IG_aflex_v2.mat']);
v3 = load([subj,'_T',num2str(trial),'_Opt7_IG_aflex_v3.mat']);

figure()
subplot(221)
plot(s.R.exp.qspline,'k','LineWidth',1.5); hold on
plot(s.R.x,'r','LineWidth',1.5); hold on
plot(v1.R.x,'Color',[168 218 220]./255,'LineWidth',1.5); hold on
plot(v2.R.x,'Color',[69 123 157]./255,'LineWidth',1.5); hold on
plot(v3.R.x,'Color',[0 180 216]./255,'LineWidth',1.5); hold on
title('T1'); box off

trial = 2;
s  = load([subj,'_T',num2str(trial),'_Opt7_IG_standard.mat']);
v1 = load([subj,'_T',num2str(trial),'_Opt7_IG_aflex_v1.mat']);
v2 = load([subj,'_T',num2str(trial),'_Opt7_IG_aflex_v2.mat']);
v3 = load([subj,'_T',num2str(trial),'_Opt7_IG_aflex_v3.mat']);

subplot(222)
plot(s.R.exp.qspline,'k','LineWidth',1.5); hold on
plot(s.R.x,'r','LineWidth',1.5); hold on
plot(v1.R.x,'Color',[168 218 220]./255,'LineWidth',1.5); hold on
plot(v2.R.x,'Color',[69 123 157]./255,'LineWidth',1.5); hold on
plot(v3.R.x,'Color',[0 180 216]./255,'LineWidth',1.5); hold on
title('T2'); box off

suptitle('TD5- IG aflex')

%% CP 2
subj  = 'CP2';
trial = 4;
s  = load([subj,'_T',num2str(trial),'_Opt7_IG_standard.mat']);
v1 = load([subj,'_T',num2str(trial),'_Opt7_IG_aflex_v1.mat']);
v2 = load([subj,'_T',num2str(trial),'_Opt7_IG_aflex_v2.mat']);
v3 = load([subj,'_T',num2str(trial),'_Opt7_IG_aflex_v3.mat']);

figure()
subplot(221)
plot(s.R.exp.qspline,'k','LineWidth',1.5); hold on
plot(s.R.x,'r','LineWidth',1.5); hold on
plot(v1.R.x,'Color',[168 218 220]./255,'LineWidth',1.5); hold on
plot(v2.R.x,'Color',[69 123 157]./255,'LineWidth',1.5); hold on
plot(v3.R.x,'Color',[0 180 216]./255,'LineWidth',1.5); hold on
title('T4'); box off

trial = 5;
s  = load([subj,'_T',num2str(trial),'_Opt7_IG_standard.mat']);
v1 = load([subj,'_T',num2str(trial),'_Opt7_IG_aflex_v1.mat']);
v2 = load([subj,'_T',num2str(trial),'_Opt7_IG_aflex_v2.mat']);
v3 = load([subj,'_T',num2str(trial),'_Opt7_IG_aflex_v3.mat']);

subplot(222)
plot(s.R.exp.qspline,'k','LineWidth',1.5); hold on
plot(s.R.x,'r','LineWidth',1.5); hold on
plot(v1.R.x,'Color',[168 218 220]./255,'LineWidth',1.5); hold on
plot(v2.R.x,'Color',[69 123 157]./255,'LineWidth',1.5); hold on
plot(v3.R.x,'Color',[0 180 216]./255,'LineWidth',1.5); hold on
title('T5'); box off

trial = 6;
s  = load([subj,'_T',num2str(trial),'_Opt7_IG_standard.mat']);
v1 = load([subj,'_T',num2str(trial),'_Opt7_IG_aflex_v1.mat']);
v2 = load([subj,'_T',num2str(trial),'_Opt7_IG_aflex_v2.mat']);
v3 = load([subj,'_T',num2str(trial),'_Opt7_IG_aflex_v3.mat']);

subplot(223)
plot(s.R.exp.qspline,'k','LineWidth',1.5); hold on
plot(s.R.x,'r','LineWidth',1.5); hold on
plot(v1.R.x,'Color',[168 218 220]./255,'LineWidth',1.5); hold on
plot(v2.R.x,'Color',[69 123 157]./255,'LineWidth',1.5); hold on
plot(v3.R.x,'Color',[0 180 216]./255,'LineWidth',1.5); hold on
title('T6'); box off

trial = 7;
s  = load([subj,'_T',num2str(trial),'_Opt7_IG_standard.mat']);
v1 = load([subj,'_T',num2str(trial),'_Opt7_IG_aflex_v1.mat']);
v2 = load([subj,'_T',num2str(trial),'_Opt7_IG_aflex_v2.mat']);
v3 = load([subj,'_T',num2str(trial),'_Opt7_IG_aflex_v3.mat']);

subplot(224)
plot(s.R.exp.qspline,'k','LineWidth',1.5); hold on
plot(s.R.x,'r','LineWidth',1.5); hold on
plot(v1.R.x,'Color',[168 218 220]./255,'LineWidth',1.5); hold on
plot(v2.R.x,'Color',[69 123 157]./255,'LineWidth',1.5); hold on
plot(v3.R.x,'Color',[0 180 216]./255,'LineWidth',1.5); hold on
title('T7'); box off

suptitle('CP2- IG aflex')

%%
subj  = 'CP4';
trial = 3;
s  = load([subj,'_T',num2str(trial),'_Opt7_IG_standard.mat']);
v1 = load([subj,'_T',num2str(trial),'_Opt7_IG_aflex_v1.mat']);
v2 = load([subj,'_T',num2str(trial),'_Opt7_IG_aflex_v2.mat']);
v3 = load([subj,'_T',num2str(trial),'_Opt7_IG_aflex_v3.mat']);

figure()
subplot(221)
plot(s.R.exp.qspline,'k','LineWidth',1.5); hold on
plot(s.R.x,'r','LineWidth',1.5); hold on
plot(v1.R.x,'Color',[168 218 220]./255,'LineWidth',1.5); hold on
plot(v2.R.x,'Color',[69 123 157]./255,'LineWidth',1.5); hold on
plot(v3.R.x,'Color',[0 180 216]./255,'LineWidth',1.5); hold on
title('T3'); box off

trial = 4;
s  = load([subj,'_T',num2str(trial),'_Opt7_IG_standard.mat']);
v1 = load([subj,'_T',num2str(trial),'_Opt7_IG_aflex_v1.mat']);
v2 = load([subj,'_T',num2str(trial),'_Opt7_IG_aflex_v2.mat']);
v3 = load([subj,'_T',num2str(trial),'_Opt7_IG_aflex_v3.mat']);

subplot(222)
plot(s.R.exp.qspline,'k','LineWidth',1.5); hold on
plot(s.R.x,'r','LineWidth',1.5); hold on
plot(v1.R.x,'Color',[168 218 220]./255,'LineWidth',1.5); hold on
plot(v2.R.x,'Color',[69 123 157]./255,'LineWidth',1.5); hold on
plot(v3.R.x,'Color',[0 180 216]./255,'LineWidth',1.5); hold on
title('T4'); box off

suptitle('CP4- IG aflex')

%% CP 8
subj  = 'CP8';
trial = 2;
s  = load([subj,'_T',num2str(trial),'_Opt7_IG_standard.mat']);
v1 = load([subj,'_T',num2str(trial),'_Opt7_IG_aflex_v1.mat']);
v2 = load([subj,'_T',num2str(trial),'_Opt7_IG_aflex_v2.mat']);
v3 = load([subj,'_T',num2str(trial),'_Opt7_IG_aflex_v3.mat']);

figure()
subplot(221)
plot(s.R.exp.qspline,'k','LineWidth',1.5); hold on
plot(s.R.x,'r','LineWidth',1.5); hold on
plot(v1.R.x,'Color',[168 218 220]./255,'LineWidth',1.5); hold on
plot(v2.R.x,'Color',[69 123 157]./255,'LineWidth',1.5); hold on
plot(v3.R.x,'Color',[0 180 216]./255,'LineWidth',1.5); hold on
title('T2'); box off

trial = 3;
s  = load([subj,'_T',num2str(trial),'_Opt7_IG_standard.mat']);
v1 = load([subj,'_T',num2str(trial),'_Opt7_IG_aflex_v1.mat']);
v2 = load([subj,'_T',num2str(trial),'_Opt7_IG_aflex_v2.mat']);
v3 = load([subj,'_T',num2str(trial),'_Opt7_IG_aflex_v3.mat']);

subplot(222)
plot(s.R.exp.qspline,'k','LineWidth',1.5); hold on
plot(s.R.x,'r','LineWidth',1.5); hold on
plot(v1.R.x,'Color',[168 218 220]./255,'LineWidth',1.5); hold on
plot(v2.R.x,'Color',[69 123 157]./255,'LineWidth',1.5); hold on
plot(v3.R.x,'Color',[0 180 216]./255,'LineWidth',1.5); hold on
title('T3'); box off

trial = 4;
s  = load([subj,'_T',num2str(trial),'_Opt7_IG_standard.mat']);
v1 = load([subj,'_T',num2str(trial),'_Opt7_IG_aflex_v1.mat']);
v2 = load([subj,'_T',num2str(trial),'_Opt7_IG_aflex_v2.mat']);
v3 = load([subj,'_T',num2str(trial),'_Opt7_IG_aflex_v3.mat']);

subplot(223)
plot(s.R.exp.qspline,'k','LineWidth',1.5); hold on
plot(s.R.x,'r','LineWidth',1.5); hold on
plot(v1.R.x,'Color',[168 218 220]./255,'LineWidth',1.5); hold on
plot(v2.R.x,'Color',[69 123 157]./255,'LineWidth',1.5); hold on
plot(v3.R.x,'Color',[0 180 216]./255,'LineWidth',1.5); hold on
title('T4'); box off

trial = 5;
s  = load([subj,'_T',num2str(trial),'_Opt7_IG_standard.mat']);
v1 = load([subj,'_T',num2str(trial),'_Opt7_IG_aflex_v1.mat']);
v2 = load([subj,'_T',num2str(trial),'_Opt7_IG_aflex_v2.mat']);
v3 = load([subj,'_T',num2str(trial),'_Opt7_IG_aflex_v3.mat']);

subplot(224)
plot(s.R.exp.qspline,'k','LineWidth',1.5); hold on
plot(s.R.x,'r','LineWidth',1.5); hold on
plot(v1.R.x,'Color',[168 218 220]./255,'LineWidth',1.5); hold on
plot(v2.R.x,'Color',[69 123 157]./255,'LineWidth',1.5); hold on
plot(v3.R.x,'Color',[0 180 216]./255,'LineWidth',1.5); hold on
title('T5'); box off

suptitle('CP8- IG aflex')

%%
subj  = 'CP10';
trial = 13;
s  = load([subj,'_T',num2str(trial),'_Opt7_IG_standard.mat']);
v1 = load([subj,'_T',num2str(trial),'_Opt7_IG_aflex_v1.mat']);
v2 = load([subj,'_T',num2str(trial),'_Opt7_IG_aflex_v2.mat']);
v3 = load([subj,'_T',num2str(trial),'_Opt7_IG_aflex_v3.mat']);

figure()
subplot(221)
plot(s.R.exp.qspline,'k','LineWidth',1.5); hold on
plot(s.R.x,'r','LineWidth',1.5); hold on
plot(v1.R.x,'Color',[168 218 220]./255,'LineWidth',1.5); hold on
plot(v2.R.x,'Color',[69 123 157]./255,'LineWidth',1.5); hold on
plot(v3.R.x,'Color',[0 180 216]./255,'LineWidth',1.5); hold on
title('T13'); box off

trial = 14;
s  = load([subj,'_T',num2str(trial),'_Opt7_IG_standard.mat']);
v1 = load([subj,'_T',num2str(trial),'_Opt7_IG_aflex_v1.mat']);
v2 = load([subj,'_T',num2str(trial),'_Opt7_IG_aflex_v2.mat']);
v3 = load([subj,'_T',num2str(trial),'_Opt7_IG_aflex_v3.mat']);

subplot(222)
plot(s.R.exp.qspline,'k','LineWidth',1.5); hold on
plot(s.R.x,'r','LineWidth',1.5); hold on
plot(v1.R.x,'Color',[168 218 220]./255,'LineWidth',1.5); hold on
plot(v2.R.x,'Color',[69 123 157]./255,'LineWidth',1.5); hold on
plot(v3.R.x,'Color',[0 180 216]./255,'LineWidth',1.5); hold on
title('T14'); box off

suptitle('CP10- IG aflex')