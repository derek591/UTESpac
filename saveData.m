function output = saveData(info,output,dataInfo,headers,tableNames, rawFlux, template)
% saveData saves output data in the site folder withinn the root folder
display('saving data');
% include dataInfo and headers in output folder
output.dataInfo = dataInfo;
output.tableNames = tableNames;

% sort structure
output = orderfields(output);

% find correct filename and save output
if strcmp(info.PF.globalCalculation,'global')
    PFtype = 'GPF_';
else
    PFtype = 'LPF_';
end
if strcmp(info.detrendingFormat,'linear')
    detrendType = 'LinDet_';
else
    detrendType = 'ConstDet_';
end
%% save .mat structure
fileName = strcat(info.siteFolder(5:end),'_',num2str(info.avgPer),'minAvg_',PFtype,detrendType,info.date,'.mat');
save(strcat(info.rootFolder,filesep,info.siteFolder,filesep,'output',filesep,fileName),'output');
%% save .csv file
if info.saveCSV
    csvSave(template,output,info)
end
%% save netCDF  From: http://stackoverflow.com/questions/21053406/matlab-save-cell-arrays-to-netcdf-file
if info.saveNetCDF
    
    fileName = strcat(info.siteFolder(5:end),'_',num2str(info.avgPer),'minAvg_',PFtype,detrendType,info.date,'.nc');
    
    ncFullFileName = strcat(info.rootFolder,filesep,info.siteFolder,filesep,'output',filesep,fileName);
    
    delete(ncFullFileName);
    
    % condition structure: convert cell strings to char arrays and logicals to doubles
    allFields = fields(output);
    for ii = 1:numel(allFields)
        
        % turn header cell strings to char arrays
        if iscell(output.(allFields{ii})) && size(output.(allFields{ii}),1)<3 && size(output.(allFields{ii}),1)>0
            output.(allFields{ii}) = output.(allFields{ii})(1,:);
            for j = 1:size(output.(allFields{ii}),2)
                output.(allFields{ii}){1,j} =  ['(',num2str(j),')',output.(allFields{ii}){1,j},' '];
            end
            output.(allFields{ii}) = cell2mat(output.(allFields{ii}))';
            
            % turn logicals to floats
        elseif islogical(output.(allFields{ii}))
            output.(allFields{ii}) = double(output.(allFields{ii}));
            
        end
    end
    
    struct2nc(output,ncFullFileName)
    ncdisp(ncFullFileName);
end
%% save raw fluxes
if info.saveRawConditionedData
    fileName = strcat(info.siteFolder(5:end),'_raw_',PFtype,detrendType,info.date);
    save(strcat(info.rootFolder,filesep,info.siteFolder,filesep,'output',filesep,fileName),'rawFlux');
end
end
