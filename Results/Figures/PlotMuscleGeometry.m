function f = PlotMuscleGeometry(R, info)
%Plot muscle geometry 

% Info on figure
lw = 1.5; % linewidth
f  = figure('Name','Muscle Geometry');

nPhases = 3; 
if nPhases > 1
    fTabGroup = uitabgroup;
end

Cs = linspecer(3);

for i=1:nPhases
    % set the name of the tab
    if nPhases>1
        file = {'Length','Force','LMT-MA'};
        tab = uitab(fTabGroup, 'Title', file{i});
        axes('parent',tab);
    end

    if i ==1
        subplot(221)
        plot(R.C.lM_ext,'k','LineWidth',lw); hold on; 
        title('lM'); ylabel('[m]'); box off
        subplot(222)
        plot(R.lMtilda,'k','LineWidth',lw); hold on
        title('lMtilda'); box off
        subplot(223)
        plot(R.C.lT_ext,'k','LineWidth',lw); hold on;
        title('lT'); ylabel('[m]'); xlabel('frames (200Hz)'); box off
        subplot(224)
        plot(R.C.lTtilda_ext,'k','LineWidth',lw); hold on;
        title('lTtilda'); xlabel('Frames (200Hz)'); box off; 
        %suptitle([info.subj ' Trial ', num2str(info.trial)])
        
    elseif i == 2
        subplot(231)
        plot(R.C.Fce,'k','LineWidth',lw); hold on
        title('Fce'); box off; 
        subplot(232)
        plot(R.C.Fpe_ext,'k','LineWidth',lw); hold on
        title('Fpe'); box off; 
        subplot(233)
        plot(R.C.fse_ext,'k','LineWidth',lw); hold on
        title('Fse'); box off; 
        subplot(234)
        plot(R.C.FM,'k','LineWidth',lw); hold on
        title('FM'); xlabel('Frames (200Hz)'); box off;
        subplot(235)
        plot(R.C.FT_ext,'k','LineWidth',lw); hold on
        title('FT'); xlabel('Frames (200Hz)'); box off
        %suptitle([info.subj ' Trial ', num2str(info.trial)])
        
    else
        subplot(211)
        plot(R.C.lMT_ext,'k','LineWidth',lw); hold on
        title('lMT'); box off;
        subplot(212)
        plot(R.C.MA_ext,'k','LineWidth',lw); hold on
        title('MA'); box off; 
    end
        
end

