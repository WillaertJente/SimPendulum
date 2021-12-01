function [error_force] = ForceEquilibrium(Fce, Fpe, lMT, lT, lM, fse)
%Force Equilibrium

FM  = Fce + Fpe;

% Force equilibrium
cos_alpha    = (lMT(1,:)-lT)./lM;
error_force  = FM.*cos_alpha - fse; 

end

