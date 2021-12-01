function [FT, lM, lT, fse, w] = CalculateTendonForce(lMtilda,lM_projected, params_OS, lMT, shift)
%Calculate Tendon Force,muscle fiber length and tendon length

% w 
lMo    = params_OS.MT(2,:); 
alphao = params_OS.MT(4,:); 
w_ext  = lMo(1).* sin(alphao(1));
w_flex = lMo(2).* sin(alphao(2));
w      = w_ext; 


% lM (muscle fiber length) 
lM_ext  = lMtilda.* lMo(1); % voor de extensor  
lM_flex = lMtilda.* lMo(2); % voor de flexor lMtilda nog aanpassen
lM      = lM_ext; 

% lT (tendon length)
lMT_ext  = lMT(1,:);
lMT_flex = lMT(2,:);
lT_ext   = lMT_ext - lM_projected;
lT_flex  = lMT_flex - lM_projected; % voor de flexor, lM projected nog aanpassen
lT       = lT_ext; 

% lTs (Tendon slack length)
lTs    = params_OS.MT(3,:); 

% lTtilda
lTtilda_ext  = lT_ext./lTs(1); 
lTtilda_flex = lT_flex./lTs(2); 

% Fse
fse_ext      = (exp(35*(lTtilda_ext - 0.995)))/5-0.25 + shift;
fse_flex     = (exp(35*(lTtilda_flex - 0.995)))/5-0.25 + shift;
fse          = fse_ext;             % Change if we want also flexor 
% FMo
FMo          = params_OS.MT(1,:);

% Compute tendon force
FT_ext       = FMo(1).* fse_ext; 
FT_flex      = FMo(2).* fse_ext; 
FT           = FT_ext; 
end

