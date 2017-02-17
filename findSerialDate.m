function [data, dataInfo, info] = findSerialDate(data, dataInfo, info)
fprintf('\nStoring serial date number in first column of each table and deleting columns 2 - 4\n')
for i = 1:length(data)
    fprintf('Finding serial dates for table %g\n',i)
    % check to make sure table loaded
    if ~isempty(data{i})
        
        % store in col 1 of table
        data{i}(:,1) = campbellDate2SerialDate(data{i});
        
        % delete columns 2-4
        data{i}(:,2:4) = [];
        
        % store begin and end date in dateInfo
        row = find(cellfun(@isempty,dataInfo(:,i))==1,1,'first');
        
        if isempty(row)
            dataInfo{end+1,i} = strcat('beg date:',datestr(data{i}(1,1)));
            dataInfo{end+1,i} = strcat('end date:',datestr(data{i}(end,1)));
        else
            dataInfo{row,i} = strcat('beg date:',datestr(data{i}(1,1)));
            dataInfo{row+1,i} = strcat('end date:',datestr(data{i}(end,1)));
        end
        
        info.date = datestr(data{i}(1,1),'yyyy_mm_dd');
        
    end
end
end

