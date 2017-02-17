function header = importHeader(filename)
% define delimter
delimiter = ',';

% Open the text file.
fileID = fopen(filename,'r');

% find number of fields
numFields = numel(strfind(fgetl(fileID),','));

% set format spec
formatSpec = strcat(repmat('%s',1,numFields),'%[^\n\r]');

% reset indicator
frewind(fileID)

% read data
header = textscan(fileID, formatSpec,'Delimiter', delimiter);

% Close the text file.
fclose(fileID);

% remove double quotes and change cell elements to type char
for ii = 1:length(header)
        header{ii} = strrep(header{ii}{:},'"','');
end