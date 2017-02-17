function [data, dataInfo, info] = date(data, dataInfo, info)
display(sprintf('\nStoring serial date number in first column of each table and deleting columns 2 - 4'))
for i = 1:length(data)
    display(sprintf('Finding serial dates for table %g',i))
    % check to make sure table loaded
    if ~isempty(data{i})
        
        % store in col 1 of table
        data{i}(:,1) = campbellDate2SerialDate(data{i});
        
        % delete columns 2-4
        data{i}(:,2:4) = [];
        
        % store begin and end date in dateInfo
        row = find(cellfun(@isempty,dataInfo(:,i))==1,1,'first');
        
        if isempty(row)
            dataInfo{end+1,i} = strcat('beg date:',datestr(date(1)));
            dataInfo{end+1,i} = strcat('end date:',datestr(date(end)));
        else
            dataInfo{row,i} = strcat('beg date:',datestr(date(1)));
            dataInfo{row+1,i} = strcat('end date:',datestr(date(end)));
        end
        
        info.date = datestr(date(1),'yyyy_mm_dd');
        
    end
end
end

