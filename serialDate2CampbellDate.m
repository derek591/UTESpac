function campbellDateVec = serialDate2CampbellDate(tSerial)

% create date vector from filled time stamps
[yy, ~, ~, hh, mm, ss] = datevec(tSerial);

% round seconds to 2 decimals to eliminate round off error
ss = round(ss*100)./100;

% find campbell style date vector and store in col 1:4
campbellDateVec = [yy', floor(tSerial') - datenum(yy,1,0)', hh'*100 + mm', ss'];
