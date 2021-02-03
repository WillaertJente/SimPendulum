function [offset] = CalculateOffset(IKtrial, params)
%Calculate offset voor bepalen LMT en MA
% x (BK) + offset = knee angle die nodig is voor lMT en ma te bepalen

IK         = importdata(IKtrial); 

col_pelT   = find(strcmp(IK.colheaders,'pelvis_tilt'));

if params.z == 18
    col_hipF   = find(strcmp(IK.colheaders,'hip_flexion_l'));
else
    col_hipF   = find(strcmp(IK.colheaders,'hip_flexion_r'));
end

start      = ones(length(IK.data),1)*90;
offset     = start - IK.data(:,col_pelT) - IK.data(:,col_hipF);
end

