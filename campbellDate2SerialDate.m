function tSerial = campbellDate2SerialDate(campbellDateVec)

tSerial = datenum(campbellDateVec(:,1),1,0) + campbellDateVec(:,2) + floor(campbellDateVec(:,3)./100)./24 + ...
    mod(campbellDateVec(:,3),floor(campbellDateVec(:,3)./100)*100)./(60*24) + campbellDateVec(:,4)./(60*60*24);