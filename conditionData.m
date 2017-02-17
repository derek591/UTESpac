function [data, output] = conditionData(data,info,tableNames,template,headers)


display(sprintf('\nConditioning Data'))

% find number of tables
numTables = numel(data);

% iterate through data tables
for i = 1:numTables
    
    if isempty(data{i}); continue; end;
    
    % initialize spikeDef and minMax arrays for spike and absolute limits tests
    tableColumns = size(data{i},2);
    spikeDef = nan(1,tableColumns);
    minMax = nan(2,tableColumns);
    
    % find name of all variables for the spike and abolute limits test
    spikeDefFields = fields(info.spikeTest.spikeDef);
    
    % iterate through all spike Fields to build definitions in spikeDef and minMax
    for j = 1:length(spikeDefFields)
        jthField = spikeDefFields{j};
        
        if strcmp(jthField,'otherInstrument'); continue; end
        
        % find columns that fall into jth field
        sensorFlag = strncmp(template.(jthField),headers{i}(1,:),length(template.(jthField))-1);
        
        % for T_Sonic, ensure that everything is in C, not K.
        if strcmp(jthField,'Tson') 
             for k = find(sensorFlag)
                 if nanmedian(data{i}(:,k)) > 250
                     data{i}(:,k) = data{i}(:,k) - 273.15;
                 end
             end
        end
        
        % update spike definition and minMax arrays
        spikeDef(sensorFlag) = info.spikeTest.spikeDef.(spikeDefFields{j});
        minMax(:,sensorFlag) = repmat(info.absoluteLimitsTest.(spikeDefFields{j})',1,sum(sensorFlag));
    end
    
    % nan values that fall above or below the min and max set values
    for j = 1:size(data{i},2)
        if ~sum(isnan(minMax(:,j)))
            data{i}(data{i}(:,j) < minMax(1,j) | data{i}(:,j) > minMax(2,j),j) = nan;
        end
    end
    
    % set 'otherInstruments'
    spikeDef(isnan(spikeDef)) = info.spikeTest.spikeDef.otherInstrument;
    
    % re-nan timestamp and sonic diagnostic
    spikeDef(1) = nan; spikeDef(strncmp(template.sonDiagnostic,headers{i}(1,:),length(template.sonDiagnostic)-1)) = nan;
    
    % 1. SPIKE REMOVAL AND SPIKE TEST
    
    % find break points
    numDays = data{i}(end,1) - data{i}(1,1);
    numBins = round(numDays)/(info.avgPer*info.spikeTest.windowSizeFraction/1440)*2+1;
    breakPoints = round(linspace(0,size(data{i},1),numBins));
    
    % create flag to indicate no remaining spikes within bin
    noSpikeFlag = false(length(breakPoints)-1,1);
    
    % create spike flag
    spikeFlag = false(size(data{i}));
    
    % create nan flag
    nanFlag = isnan(data{i});
    
    % iterate through max number of runs
    for j = 1:info.spikeTest.maxRuns
        
        % if no spikes remain anywhere, continue to next table
        if min(noSpikeFlag(1:end-1)); continue; end;
        
        tic
        if j==1; fprintf('\nDespiking table %g - %s',i,tableNames{i}); end
        
        % iterate through date bins
        for k = 1:length(breakPoints)-2
            
            % skip bin if no spikes were found on the previous run within the bin
            if noSpikeFlag(k); continue; end
            
            % find starting row
            startRow = breakPoints(k)+1;
            
            % find end row, +2 because we iterate at 1/2 at a time (e.g. breakpoints are doubled)
            endRow = breakPoints(k+2);
            
            % store local data
            localData = data{i}(startRow:endRow,:);
            
            % store local nan flag
            localNanFlag = nanFlag(startRow:endRow,:);
            
            % find the standard deviation for each column
            stds = nanstd(localData);
            
            % find the deviations from the mean for all locl data
            deviations = bsxfun(@minus,localData,nanmean(localData));
            
            % divide all deviations by the respective standard deviation
            normalizedDeviations = abs(bsxfun(@rdivide,deviations,stds));
            
            % flag all data that falls too many standard deviations from the mean
            potSpikeFlags = bsxfun(@ge,normalizedDeviations,spikeDef);
            
            % remove consecutive spike flags that are considered to be possibly physical
            localSpikeFlags = consecFlagRemoval(potSpikeFlags,info.spikeTest.maxConsecutiveOutliers);
            
            % store local run flags in global spike flag
            spikeFlag(startRow:endRow,:) = spikeFlag(startRow:endRow,:) | localSpikeFlags;
            
            % if there are no localSpikeFlags, change the kth noSpikeFlag and continue
            if ~any(any(localSpikeFlags))
                noSpikeFlag(k) = 1;
                continue
            end
            
            % --- linearly interpolate over spike values in localData
            for col=2:tableColumns
                
                % flag all spiked and nan'd values
                badFlag = localSpikeFlags(:,col) | localNanFlag(:,col);
                
                % ensure that unspiked values exist before interpolating, otherwise nan full period
                if sum(badFlag)/size(badFlag,1) > 0.6
                    localData(:,col) = localData(:,col)*nan;
                    continue
                end
                
                % create array with dx = 1
                X = 1:size(badFlag,1);
                % store local column in Y
                Y = localData(:,col);
                
                % queried values are those marked as a spike
                Xq = X(localSpikeFlags(:,col));
                
                % replace spiked values with interp1
                Y(Xq) = interp1(X(~badFlag),Y(~badFlag),Xq);
                
                localData(:,col) = Y;
            end
            
            % put local data back into data cell
            data{i}(startRow:endRow,:) = localData;
        end
        fprintf('Pass %g time: %g seconds\n',j,toc)
    end
    
    %--- create mean spike flags and nan flags
    meanSpikeFlag = simpleAvg([data{i}(:,1), spikeFlag(:,2:end)],info.avgPer); % find average num of spike flags per avg period
    meanSpikeFlag(:,1) = 0; % set date column to 0
    meanSpikeFlag(meanSpikeFlag < info.spikeTest.maxPercent/100) = 0;  % zero percentages less than maxPercent
    output.([tableNames{i},'SpikeFlag']) = logical(meanSpikeFlag); % create logical with remainder
    
    meanNanFlag = simpleAvg([data{i}(:,1), nanFlag(:,2:end)],info.avgPer); % find average num of Nan flags per avg period
    meanNanFlag(:,1) = 0; % set date column to 0
    meanNanFlag(meanNanFlag < info.nanTest.maxPercent/100) = 0;  % zero percentages less than maxPercent
    output.([tableNames{i},'NanFlag']) = logical(meanNanFlag); % create logical with remainder
end