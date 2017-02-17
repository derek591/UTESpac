function [headersCell, dataFiles, tableNames, info] = findFiles(info)

% find possible sites of interest
sitesStruct = dir(strcat(info.rootFolder,filesep,'site*'));

% display possible sites to command window
for ii = 1:numel(sitesStruct)
    display(sprintf('%g. %s',ii,sitesStruct(ii).name));
end

% ask user to select appropriate site
site = sitesStruct(input('Please indicate site number of interest: ')).name; clc;
info.siteFolder = site;

% store site folder
siteFolder = strcat(info.rootFolder,filesep,site);

% run siteInfo script
run(strcat(siteFolder,filesep,'siteInfo.m'));

% find headers for selected site
headers = dir(strcat(siteFolder,filesep,'*header.*'));

if isempty(headers)
    error('Unable to Find Table Headers!  Check format of header name.')
end

% iterate through headers and store .csv files in csvFileStruct
for ii = 1:length(headers)
    currentHeaderFileName = headers(ii).name;
    try
        % find current header array
        currentHeaderArray =  importHeader(strcat(info.rootFolder,filesep,site,filesep,currentHeaderFileName));
               
        % store measurement height on second row of currentHeaderArray
        for jj = 1:length(currentHeaderArray)
            
            % find digits in column title
            digitLocations = find(isstrprop(currentHeaderArray{1,jj},'digit')==1);
            
            % find decimals in column title
            decimalLocations = strfind(currentHeaderArray{1,jj},'.');
            
            % concatnate digits and decimals, sort in ascending order
            digitAndDecimalLocations = sort(horzcat(digitLocations, decimalLocations));
            
            % if digitsAndDecimals is empty, no height was found 
            if isempty(digitAndDecimalLocations)
                currentHeaderArray{2,jj} = [];
                
                % if digitsAndDecimals has 1 value, this must be the height 
            elseif numel(digitAndDecimalLocations) == 1
                heightLocations = digitAndDecimalLocations;
                currentHeaderArray{2,jj} = str2double(currentHeaderArray{1,jj}(heightLocations));
                
            else % if digitsAndDecimals has multiple values, use diff to find consecutive locations, 
                 % last chunk of consecutive locations must be the height locations
                heightBeginMarker = find(diff(digitAndDecimalLocations) > 1,1,'last');
                if isempty(heightBeginMarker);
                    heightBeginMarker = 0;
                end
                heightLocations = digitAndDecimalLocations(heightBeginMarker+1:end);
                currentHeaderArray{2,jj} = str2double(currentHeaderArray{1,jj}(heightLocations));
            end
        end
        
        
        % store headers in cell 
        headersCell{ii} = currentHeaderArray; clear currentHeaderArray
        
    catch err
        error('Unable to load %s \@ ln %g.  Check format! \n%s',currentHeaderFileName,err.message);
    end
    
    % find length of table name
    tableNameLength = strfind(currentHeaderFileName,'_header')-1;
    
    % remove '_header' to find table name
    tableName = currentHeaderFileName(1:tableNameLength);
    
    % store table names in cell
    tableNames{ii} = tableName;
    
    % find all files associated with the tableName
    tableFiles = dir(strcat(siteFolder,filesep,strcat('*',tableName,'*')));
    
    % eliminate the header file from the list
    tableFiles(~cellfun(@isempty,strfind({tableFiles(:).name},'header'))) = [];
    
    % sort and store csv files in csvFilesCell
    csvFiles = cell(length(tableFiles),1);  %preallocate for speed
    csvDates = NaN(size(csvFiles));
    for jj = 1:length(tableFiles)
        tableFile = tableFiles(jj).name;
        
        % place tableFile name in csvFiles cell
        csvFiles{jj} = tableFile;
        
        % find location of date
        dateBeginLoc = strfind(tableFiles(jj).name,tableName)+length(tableName);

        % parse date string to create date number
        csvDates(jj) = datenum(str2double(tableFile(dateBeginLoc+1:dateBeginLoc+4)),str2double(tableFile(dateBeginLoc+6:dateBeginLoc+7)),... % year, month
            str2double(tableFile(dateBeginLoc+9:dateBeginLoc+10)),str2double(tableFile(dateBeginLoc+12:dateBeginLoc+13)),... % day, hour
            str2double(tableFile(dateBeginLoc+14:dateBeginLoc+15)),0); % minute, second
    end
    
    % place csv files in cell matrix where the row is determined by the date
    if ii == 1
        dateBegin = floor(min(csvDates)-20); % all tables must begin and end within 20 days of first and last dates of table1
        dateEnd = floor(max(csvDates)+20);
        dataFiles = cell(dateEnd-dateBegin,length(headers));
    end
    for jj = 1:length(csvFiles)
        dataFiles{floor(csvDates(jj))-dateBegin,ii} = csvFiles{jj};
    end
end
% check cell matrix for empty rows and delete them
%dataFiles(cellfun(@isempty,dataFiles(:,1)) & cellfun(@isempty,dataFiles(:,2)),:) = [];
dataFiles(logical(sum(cellfun(@isempty,dataFiles),2))) = [];

% ask user to select which dates to calculate
displayCell = cell(size(dataFiles,1),size(dataFiles,2)+1);
displayCell(1:end,1) = num2cell(1:size(dataFiles,1))';
displayCell(1:end,2:end) = dataFiles;

display(sprintf('\nDisplaying all files'))
showcell(displayCell);

SOI = input('Plese input dates of interest. e.g. [1 3 4:7] or ''0'' for all dates: ');

if SOI > 0
    dataFiles = dataFiles(SOI,:);
end
for ii = 1:size(dataFiles,2)
    for jj = 1:size(dataFiles,1)
        dataFiles{jj,ii} = strcat(info.rootFolder,filesep,site,filesep,dataFiles{jj,ii});
    end
end
clc
end


