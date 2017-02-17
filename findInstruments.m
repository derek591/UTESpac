function [sensorInfo, info] = findInstruments(headers,template,info)

% extract variable names from templateStruct
templateVariables = fields(template);

% iterate through headers
for ii = 1:numel(headers)
    
    % find current header
    currentHeader = headers{ii}(1,:);
    
    % loop through all variable names for template matches
    for jj = 1:length(templateVariables)
        
        % find current template in the loop
        currentTemplate = template.(templateVariables{jj});
        
        % check to make sure the template is not empty
        if isempty(currentTemplate)
            continue
        end
        
        % find matches for current header and template
        currentMatches = strfndw(currentHeader,currentTemplate);
        
        % find number of current matches
        numMatches = numel(currentMatches);
        
        % find match heights
        currentMatchHeights = nan(size(currentMatches));
        for kk = 1:numel(currentMatches)
            tempMatchHeight = sscanf(currentHeader{currentMatches(kk)}(strfind(currentTemplate,'*'):end),'%f');
            if ~isempty(tempMatchHeight)
                currentMatchHeights(kk) =   tempMatchHeight;
            end
        end
        
        % if matches exist store them, along with height, in sensorInfo.  Column 1 is table number, Column 2 is
        % column number with in the table and column 3 is the height.
        if ~isempty(currentMatches)
            if ~exist('sensorInfo','var')
                sensorInfo = struct;
            end
            if isfield(sensorInfo,sprintf('%s',templateVariables{jj}))
                sensorInfo.(templateVariables{jj})(end+1:end+numMatches,1) = ii;
                sensorInfo.(templateVariables{jj})(end-numMatches+1:end,2) = currentMatches;
                sensorInfo.(templateVariables{jj})(end-numMatches+1:end,3) = currentMatchHeights;
            else
                sensorInfo.(templateVariables{jj})(1:numMatches,1) = ii;
                sensorInfo.(templateVariables{jj})(1:numMatches,2) = currentMatches;
                sensorInfo.(templateVariables{jj})(1:numMatches,3) = currentMatchHeights;
            end
        end
    end
end

if ~exist('sensorInfo','var')
    sensorInfo = struct();
end
% output number of found instruments
for ii = 1:length(templateVariables)
    variable = templateVariables{ii};
    if isfield(sensorInfo,variable)
        numSensors = size(sensorInfo.(templateVariables{ii}),1);
    else
        numSensors = 0;
    end
    display(sprintf('variable name: %s.  Sensors found: %g',variable,numSensors));
end

flag1 = input('\nIs the number of discovered sensors correct(1 for yes, 0 for no): ');
if ~flag1
    error('Check templates!')
end
clc
% output sonic details
if isfield(sensorInfo,'u')
    numSonics = size(sensorInfo.u,1);
    for ii = 1:numSonics
        height = sensorInfo.u(ii,3);
        orientation = info.sonicOrientation(ii);
        if info.sonicManufact(ii)
            manufact = 'Campbell';
        else
            manufact = 'RMYoung ';
        end
        sensorInfo.u(ii,4) = orientation;
        sensorInfo.u(ii,5) = info.sonicManufact(ii);
        display(sprintf('Sonic%g: Manufacturer=%s -  height=%gm - Orientation=%gdeg',ii,manufact,height,orientation))
    end
    flag2 = input('\nIs the Sonic Information Correct (1 for yes, 0 for no): ');
    if ~ flag2
        error('Program stopped by user.  Check header format and height information')
    end
    % store sonic height in info struct
    info.sonicHeight = sensorInfo.u(:,3)';
else
    flag3 = input('\nNo sonics found.  Is this correct? (1 for yes, 0 for no): ');
    if ~ flag3
        error('Program stopped by user.  Check header format and height information')
    end
end

info = orderfields(info);
clc
end

