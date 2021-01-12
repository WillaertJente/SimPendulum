function [offset] = CalculateOffset(IKtrial)
%Calculate offset voor bepalen LMT en MA
% x (BK) + offset = knee angle die nodig is voor lMT en ma te bepalen

IK         = importdata(IKtrial); 

col_pelT   = find(strcmp(IK.colheaders,'pelvis_tilt'));
col_hipF   = find(strcmp(IK.colheaders,'hip_flexion_l'));

start      = ones(length(IK.data),1)*90;
offset     = start - IK.data(:,col_pelT) - IK.data(:,col_hipF);
end

