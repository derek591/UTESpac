function output = avg(data,info,tableNames,output,headers,sensorInfo)

for ii = 1:length(data)
    fprintf('\naveraging %s',tableNames{ii})
    currentTable = data{ii};
    if ~isempty(currentTable)
        
        % find table header
        header = headers{ii};
        
        % find table frequency
        tableFreq = 1/((currentTable(2,1) - currentTable(1,1))*(3600*24));
        
        % find number of samples per averaging period
        numSamplesPerPeriod = round(info.avgPer*60*tableFreq);
        
        % find number of averaging periods
        numAvgPeriods = size(currentTable,1)/numSamplesPerPeriod;
        
        % preallocate avgTable
        avgTable = nan(numAvgPeriods,size(currentTable,2));
        
        % iterate through avgTable rows
        for jj = 1:size(avgTable,1)
            
            % store local data
            localData = currentTable((jj-1)*numSamplesPerPeriod+1:jj*numSamplesPerPeriod,:);
            
            % average current data and store in avgTable
            avgTable(jj,:) = nanmean(localData);
            
            % date is last value of averaging period
            avgTable(jj,1) = localData(end,1);  
            
            % check for wind direction from wind bird and find mean unit vector direction
            if isfield(sensorInfo,'birdDir')
                for k = 1:length(sensorInfo.birdDir(:,1))
                    if sensorInfo.birdDir(k,1) == ii
                        WS = localData(:,sensorInfo.birdSpd(k,2));
                        WD = localData(:,sensorInfo.birdDir(k,2)).*pi./180;
                        v = nanmean(WS.*sin(WD));
                        u = nanmean(WS.*cos(WD));
                        avgTable(jj,sensorInfo.birdDir(k,2)) = mod(atan2(v,u)*180/pi,360);
                    end
                end
            end
            
        end
        tableName = tableNames{ii};
        output.(tableName) = avgTable;
        output.([tableName,'Header']) = header;
    end
end
end