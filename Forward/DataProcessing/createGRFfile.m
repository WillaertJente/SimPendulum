function [grfFileOut] = createGRFfile(fileLocation,trialName,subjName)
% Create GRF input file for inverse dynamics in Opensim


% load([fileLocation 'experimentalData\' trialName '.mat'])
load([fileLocation trialName '.mat'])

if exist('data')
    
    % put grf and cop into opensim/simm coordinate system
    grf_r = [data.analog.plateforces(:,2), -data.analog.plateforces(:,3), data.analog.plateforces(:,1)];
    grf_l = [data.analog.plateforces(:,8), -data.analog.plateforces(:,9), data.analog.plateforces(:,7)];
    
    cop_r_temp1 = [data.calculated.cop(:,4), zeros(length(data.calculated.cop),1),data.calculated.cop(:,3)];
    cop_l_temp1 = [data.calculated.cop(:,6), zeros(length(data.calculated.cop),1),data.calculated.cop(:,5)];
    
    % convert from atime to mtime
    for i=1:3
        grf_r_out(:,i) = spline(atime',grf_r(:,i),mtime');
        grf_l_out(:,i) = spline(atime',grf_l(:,i),mtime');
        cop_r_temp(:,i) = spline(atime',cop_r_temp1(:,i),mtime');
        cop_l_temp(:,i) = spline(atime',cop_l_temp1(:,i),mtime');
        if i<3
            LVDT_temp(:,i) = spline(atime',LVDT(:,i),mtime');
        end
    end
    
    mom_r_out = zeros(length(grf_r_out),3);
    mom_l_out = zeros(length(grf_r_out),3);  
  
    
else
    
    % put grf and cop into opensim/simm coordinate system
    grf_r = [GRF(:,2), GRF(:,3), GRF(:,1)];
    grf_l = [GRF(:,8), GRF(:,9), GRF(:,7)];
    
    cop_r_temp1 = [COP(:,3), zeros(length(COP),1),COP(:,4)];
    cop_l_temp1 = [COP(:,1), zeros(length(COP),1),COP(:,2)];
    
    % convert from atime to mtime
    for i=1:3
        grf_r_out(:,i) = spline(atime',grf_r(:,i),mtime');
        grf_l_out(:,i) = spline(atime',grf_l(:,i),mtime');
        cop_r_temp(:,i) = spline(atime',cop_r_temp1(:,i),mtime');
        cop_l_temp(:,i) = spline(atime',cop_l_temp1(:,i),mtime');
        if i<3
            LVDT_temp(:,i) = spline(atime',LVDT(:,i),mtime');
        end
    end
    
    mom_r_out = zeros(length(grf_r_out),3);
    mom_l_out = zeros(length(grf_r_out),3);
    
    
end


% convert cop wrt to heel marker
cop_r = cop_r_temp;
cop_l = cop_l_temp;


cop_r2(:,1) = cop_r(:,1) + LVDT_temp(:,2)*10;
cop_l2(:,1) = cop_l(:,1) + LVDT_temp(:,2)*10;

cop_r2(:,3) = cop_r(:,3) + LVDT_temp(:,1)*10;
cop_l2(:,3) = cop_l(:,3) + LVDT_temp(:,1)*10;



% move in to global frame

%%% Geometry of forceplates [Dimensions in mm]
%%% (looking from the control station)
% Global origin to Right force plate (looking from the control station)
origin_right_x = 464+464/2;
origin_right_y = 507/2;
% Global origin to Left force plate origin
origin_left_x = 464/2;
origin_left_y = 507/2;


cop_r(:,1) = cop_r2(:,1) + origin_right_y;
cop_l(:,1) = cop_l2(:,1) + origin_left_y;

cop_r(:,3) = cop_r2(:,3) + origin_right_x;
cop_l(:,3) = cop_l2(:,3) + origin_left_x;


% conver to m
cop_r = cop_r./1000;
cop_l = cop_l./1000;


grfFileOut = [mtime',grf_r_out,cop_r,grf_l_out,cop_l,mom_r_out,mom_l_out];


%%% Write data to file %%%

%   header
[fileLocation 'dataForSimulation\' subjName '_' trialName '_GRFinputFile.mot']
fid = fopen([fileLocation 'dataForSimulation\' subjName '_' trialName '_GRFinputFile.mot'],'W');
fprintf(fid,'GRFinputFile\n');
fprintf(fid,'version=1\n');
fprintf(fid,'nRows=%i\n',size(grfFileOut,1));
fprintf(fid,'nColumns=%i\n',size(grfFileOut,2));
fprintf(fid,'inDegrees=yes\n');
fprintf(fid,'endheader\n');
fprintf(fid,'time\tground_force_vx\tground_force_vy\tground_force_vz\tground_force_px\tground_force_py\tground_force_pz\t');
fprintf(fid,'1_ground_force_vx\t1_ground_force_vy\t1_ground_force_vz\t1_ground_force_px\t1_ground_force_py\t1_ground_force_pz\t');
fprintf(fid,'ground_torque_x\tground_torque_y\tground_torque_z\t1_ground_torque_x\t1_ground_torque_y\t1_ground_torque_z\n');


for j=1:size(grfFileOut,1)
    for i=1:size(grfFileOut,2)
        fprintf(fid,'%f\t',grfFileOut(j,i));
    end
    fprintf(fid,'\n');
end

fclose(fid)

