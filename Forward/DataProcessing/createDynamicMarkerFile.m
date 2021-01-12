function createDynamicMarkerFile(fileLocation,trialName,subjName)

frameRate = 120;
dataStart = 1;
numMarkers = 25;
% markers = {'C7','RASI','LASI','SACR','RKNE','RANK','RHEE','RTOE','LKNE','LANK','LHEE','LTOE'};
markers = {'LFHD','RFHD','LBHD','RBHD','C7','CLAV','RBAK','LSHO','RSHO', ...
    'RASI','LASI','RPSI','LPSI','RTHI','RKNE','RTIB','RANK','RHEE','RTOE', ...
    'LTHI','LKNE','LTIB','LANK','LHEE','LTOE'};

% load file
load([fileLocation  trialName '.mat']);


if exist('data')
    % get number of frames
    numFrames = size(data.video.markers,1);    
    
    % calculate mass (in kg)
    totalVertGRF = data.analog.plateforces(1:numFrames*9,3) + data.analog.plateforces(1:numFrames*9,9);
    meanVertGRF = nanmean(totalVertGRF)/9.81;
    disp(['Subject mass = ' num2str(meanVertGRF) ' kg']);
    
    %  Calculate SACR marker position: based on average of right and left PSIS markers    
    markerSACR = (data.video.markers(:,strmatch('RPSI',data.video.markersid,'exact'),:) + ...
        data.video.markers(:,strmatch('LPSI',data.video.markersid,'exact'),:))/2;
    
    
    %  Write data to file
    
    %   header
    fid = fopen([fileLocation 'dataForSimulation\' subjName '_' trialName '_markerFile.trc'],'W');
    fprintf(fid,'PathFileType\t4\t(X/Y/Z)\ttest.trc\n');
    fprintf(fid,'DataRate\tCameraRate\tNumFrames\tNumMarkers\tUnits\tOrigDataRate\tOrigDataStartFrame\tOrigNumFrames\n');
    fprintf(fid,'%i\t%i\t%i\t%i\tmm\t%i\t%i\t%i\n',frameRate,frameRate,numFrames,numMarkers,frameRate,dataStart,numFrames);
    fprintf(fid,'Frame#\tTime\t');
    for i=1:numMarkers
        fprintf(fid,'%s\t\t\t',markers{i});
    end
    fprintf(fid,'\n\t\t');
    for i=1:numMarkers
        fprintf(fid,'X%i\tY%i\tZ%i\t',i,i,i);
    end
    fprintf(fid,'\n\n');
    
    
    for i=1:numFrames
        time = 0.0+(i-1)*(1/120);
        fprintf(fid,'%i\t%f\t',i,time);
        for j=1:numMarkers
            if strcmp(markers{j},'SACR')
                fprintf(fid,'%f\t%f\t%f\t',markerSACR(i,2),markerSACR(i,3),markerSACR(i,1));
            else
                fprintf(fid,'%f\t%f\t%f\t',data.video.markers(i,strmatch(markers{j},data.video.markersid,'exact'),2),...
                    data.video.markers(i,strmatch(markers{j},data.video.markersid,'exact'),3),data.video.markers(i,strmatch(markers{j},data.video.markersid,'exact'),1));
            end
        end
        fprintf(fid,'\n');
    end
    
    fclose(fid);
    
else % old data format
    
    % get number of frames
    numFrames = size(Markers,1);
    
    MarkerID = cellstr(MarkerID);
    
    
    % calculate mass (in kg)
    totalVertGRF = GRF(1:numFrames*9,3) + GRF(1:numFrames*9,9);
    meanVertGRF = nanmean(totalVertGRF)/9.81;
    disp(['Subject mass = ' num2str(meanVertGRF) ' kg']);
    
    %  Calculate SACR marker position
    %   based on average of right and left PSIS markers
    
    markerSACR = (Markers(:,strmatch('RPSI',MarkerID,'exact'),:) + ...
        Markers(:,strmatch('LPSI',MarkerID,'exact'),:))/2;
    
    
    %  Write data to file
    
    %   header
    fid = fopen([fileLocation 'dataForSimulation\' subjName '_' trialName '_markerFile.trc'],'W');
    fprintf(fid,'PathFileType\t4\t(X/Y/Z)\ttest.trc\n');
    fprintf(fid,'DataRate\tCameraRate\tNumFrames\tNumMarkers\tUnits\tOrigDataRate\tOrigDataStartFrame\tOrigNumFrames\n');
    fprintf(fid,'%i\t%i\t%i\t%i\tmm\t%i\t%i\t%i\n',frameRate,frameRate,numFrames,numMarkers,frameRate,dataStart,numFrames);
    fprintf(fid,'Frame#\tTime\t');
    for i=1:numMarkers
        fprintf(fid,'%s\t\t\t',markers{i});
    end
    fprintf(fid,'\n\t\t');
    for i=1:numMarkers
        fprintf(fid,'X%i\tY%i\tZ%i\t',i,i,i);
    end
    fprintf(fid,'\n\n');
    
    
    for i=1:numFrames
        time = 0.0+(i-1)*(1/120);
        fprintf(fid,'%i\t%f\t',i,time);
        for j=1:numMarkers
            if strcmp(markers{j},'SACR')
                fprintf(fid,'%f\t%f\t%f\t',markerSACR(i,2),markerSACR(i,3),markerSACR(i,1));
            else
                fprintf(fid,'%f\t%f\t%f\t',Markers(i,strmatch(markers{j},MarkerID,'exact'),2),...
                    Markers(i,strmatch(markers{j},MarkerID,'exact'),3),Markers(i,strmatch(markers{j},MarkerID,'exact'),1));
            end
        end
        fprintf(fid,'\n');
    end
    
    fclose(fid)
    
end