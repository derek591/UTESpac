function flag = consecFlagRemoval(flag,consecRows)

% ensure that their are more rows per avg period than the consecRows option
if consecRows > size(flag,1)
    return
end

% add buffer at begining and end to prevent looping effect
flag(end+1,:) = false; flag = [false(1,size(flag,2)); flag];

flagFlag = flag;

% remove flags that don't have consecutive occurences greater than consecRows
for i = 1:consecRows
    flagFlag = flagFlag & [flag(i+1:end,:); flag(1:i,:)];
end

% remove first and last row
flag = flag(2:end-1,:);
flagFlag = flagFlag(2:end-1,:);

% add 1's to mark all consecutiveFlags greater than consecRows
for i=1:consecRows
    flagFlag = flagFlag | [flagFlag(end,:); flagFlag(1:end-1,:)];
end

% delete all consecutive flags greater than consecRows
flag(logical(flagFlag)) = false;


