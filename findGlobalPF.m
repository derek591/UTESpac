function PFinfo = findGlobalPF(info,template,sensorInfo)
% findGlobalPF loads or evaluates global PF coefficients (b0, b1, b2) based on the existence of previously calculated
% coefficients and the info.PFrange option.  For a global calculation, 5 minute local planar fit data (suffix: LPF) must
% have already been run. 

% check to see if PF needs to be calculate.  Else, load existing
if info.PF.recalculateGlobalCoefficients || ~exist([info.rootFolder,filesep,info.siteFolder,filesep,'pfInfo.mat'],'file')
    % load data
    disp('Finding Global Planar Fit Coefficients (b0, b1, b2)')
    allSites = dir([info.rootFolder,filesep,'site*']); allSites = {allSites(:).name};
    siteNum = find(strcmp(allSites,info.siteFolder));
    data = getData(info.rootFolder,'site',siteNum,'avgPer',5,'qualifier','LPF','rows',0);
    siteName = info.siteFolder(5:end);

    % find sonic heights
    for ii = 2:length(data.spdAndDirHeader)
        z(ii) = sscanf(data.spdAndDirHeader{ii},'%f');
    end
    z = unique(z); z = z(2:end);
    
    % find raw table names
    tableNames = data.tableNames;
    
    % iterate through heights
    for ii = 1:length(z)
        
        % allow user to skip height
        skipFlag = input(sprintf('Skip %g m (1=no, 0=yes)',z(ii)));
        if ~skipFlag && exist([info.rootFolder,filesep,info.siteFolder,filesep,'pfInfo.mat'],'file')
            oldPFinfo = load([info.rootFolder,filesep,info.siteFolder,filesep,'pfInfo.mat']);
            PFinfo.(['cm_',num2str(round(z(ii)*100))]) = oldPFinfo.PFinfo.(['cm_',num2str(round(z(ii)*100))]);
            continue;
        elseif ~skipFlag 
            disp('I couldn''t find an existing pfInfo.mat file so you cannot skip any levels')
        end
        
        % find variable name at ith height
        varName_u = strrep(template.u,'*',num2str(z(ii)));
        varName_v = strrep(template.v,'*',num2str(z(ii)));
        varName_w = strrep(template.w,'*',num2str(z(ii)));
        varName_sonDiag = strrep(template.sonDiagnostic,'*',num2str(z(ii)));
        
        % find table containing variable names
        flag = [];
        for j = 1:length(tableNames)
            flag = sum(strcmp(data.(strcat(tableNames{j},'Header'))(1,:),varName_u));
            localSonicTable = tableNames{j};
            if flag
                break
            end
        end
        
        % throw error if no table found, else, find appropriate columns
        if isempty(flag)
            error('Unable to find sonic table for sonic at %g m in findGloablPF.m @ ln 30-37.',z(ii))
        end
        
        % find velocity columns within correct table
        uCol = strncmp(data.(strcat(localSonicTable,'Header'))(1,:),varName_u,numel(varName_u));
        vCol = strncmp(data.(strcat(localSonicTable,'Header'))(1,:),varName_v,numel(varName_v));
        wCol = strncmp(data.(strcat(localSonicTable,'Header'))(1,:),varName_w,numel(varName_w));
        sonDiagnosticCol = strncmp(data.(strcat(localSonicTable,'Header'))(1,:),varName_sonDiag,numel(varName_sonDiag));
        spdCol = strcmp(data.spdAndDirHeader,[num2str(z(ii)),'m speed']);
        dirFlagCol = strncmp(data.spdAndDirHeader,[num2str(z(ii)),'m flag'],length([num2str(z(ii)),'m flag']));
        
        % find flags
        uNanFlag = data.(strcat(localSonicTable,'NanFlag'))(:,uCol);
        uSpikeFlag = data.(strcat(localSonicTable,'SpikeFlag'))(:,uCol);
        uWindFlag = data.spdAndDir(:,dirFlagCol);
        uTotalFlag = uNanFlag | uSpikeFlag | uWindFlag;
        vNanFlag = data.(strcat(localSonicTable,'NanFlag'))(:,vCol);
        vSpikeFlag = data.(strcat(localSonicTable,'SpikeFlag'))(:,vCol);
        vWindFlag = data.spdAndDir(:,dirFlagCol);
        vTotalFlag = vNanFlag | vSpikeFlag | vWindFlag;
        wNanFlag = data.(strcat(localSonicTable,'NanFlag'))(:,wCol);
        wSpikeFlag = data.(strcat(localSonicTable,'SpikeFlag'))(:,wCol);
        wWindFlag = data.spdAndDir(:,dirFlagCol);
        wTotalFlag = wNanFlag | wSpikeFlag | wWindFlag;
        if sum(sonDiagnosticCol>0) % ensure a diagnostic exists for sonic, else make flag all false
            sonDiagFlag = data.(localSonicTable)(:,sonDiagnosticCol); 
            sonDiagFlag(isnan(sonDiagFlag)) = 0; % false any nans in the diagnostic
            sonDiagFlag(sonDiagFlag < info.diagnosticTest.meanSonicDiagnosticLimit) = 0; % turn diangostics below the limit to 0
            sonDiagFlag = logical(sonDiagFlag);  % turn any remaining diagnostics to flagged values
        else
            sonDiagFlag = false(size(wTotalFlag));
        end       
        totalFlag = uTotalFlag | vTotalFlag | wTotalFlag | sonDiagFlag;
        
        % store variables
        u = data.(localSonicTable)(~totalFlag,uCol);
        v = data.(localSonicTable)(~totalFlag,vCol);
        w = data.(localSonicTable)(~totalFlag,wCol);
        t = data.(localSonicTable)(~totalFlag,1);
        spd = data.spdAndDir(~totalFlag,spdCol);
        
        
        % NaN wind speeds above and below info.globalPFmaxWind and info.globalPFminWind
        u(spd>info.PF.globalCalcMaxWind | spd<info.PF.globalCalcMinWind) = nan;
        v(spd>info.PF.globalCalcMaxWind | spd<info.PF.globalCalcMinWind) = nan;
        w(spd>info.PF.globalCalcMaxWind | spd<info.PF.globalCalcMinWind) = nan;
        spd(spd>info.PF.globalCalcMaxWind | spd<info.PF.globalCalcMinWind) = nan;        
        
        % step variables down to user selected, PF averaging period
        u = stpDn(u,round(info.PF.avgPer/5));
        v = stpDn(v,round(info.PF.avgPer/5));
        w = stpDn(w,round(info.PF.avgPer/5));
        t = stpDn(t,round(info.PF.avgPer/5));
        spd = stpDn(spd,round(info.PF.avgPer/5));    

        % check sonic manufacturere
        if sensorInfo.u(ii,5)
            uDir = u;
            vDir = v;
        else % for RMYoung swap u and v and multiply v by -1
            uDir = v;
            vDir = u*-1;
        end
        bearing = sensorInfo.u(ii,4);
        
        % find local direction and speed
        direction = mod(atan2(-1*vDir,uDir)*180/pi+bearing,360);
        
        % plot seasonal wind histogram and have user select sectors
        figure('units','normalized','outerposition',[0 0 1 1])
        % histogram
        subplot(2,1,1)
        hist(direction,100)
        xlim([0 360])
        ylabel('5 min occurences')
        xlabel('direction')
        changeFig(18,2,5)
        lims = ylim;
        hold on
        % time series
        subplot(2,1,2)
        plot(t,u,'b.',t,v,'g.',t,w,'r.',t,spd,'c.','markersize',10)
        title([siteName,' ',num2str(z(ii)),' m'])
        legend('u','v','w','mag')
        datetick
        axis tight

        % allow user to select binning sectors
        selectedDirection = 1;clc;
        bins = [];
        
        while ~isempty(selectedDirection)
            selectedDirection = input('Select up to 8 bin barriers, when done input ''[]'': ');
            if ~isempty(selectedDirection)
                bins(end+1) = selectedDirection;
                subplot(2,1,1)
                hold on
                plot([selectedDirection, selectedDirection], [lims(1), lims(2)],'g--','linewidth',3)
            end
        end
        clc
        bins = sort(bins);
        bins(bins>360 | bins<0) = [];
        
        % allow user to partition planar fit date windows
        subplot(2,1,2)
        hold on
        lims = ylim;
        selectedDate = 1;clc;
        dateBarriers = [];
        
        while ~isempty(selectedDate)
            selectedDate = input('Option to split data by date for different PF calculations.  Enter ''datenum(yyyy,mm,dd)'' for barrier and ''[]'' when finished: ');
            if ~isempty(selectedDate)
                dateBarriers(end+1) = selectedDate;
                plot([selectedDate, selectedDate], [lims(1), lims(2)],'g--','linewidth',3)
            end
        end
        clc
        dateBarriers = sort(dateBarriers);
        dateBarriers(dateBarriers<t(1) | dateBarriers>t(end)) = [];
        
        % find all possible days
        allDays = unique(floor(t));
        allDays(isnan(allDays)) = [];
        
        if isempty(dateBarriers)
            dateBarriers = [0, 0];
        else
            dateBarriers = [allDays(1), dateBarriers, allDays(end)];
        end
        
        % iterate through date barriers, allow user to use all days or selectively chose them
        for j = 1:length(dateBarriers)-1
            
            % find days in envelope
            if length(dateBarriers) == 2
                localDays = allDays;
            else
                localDays = allDays(find(dateBarriers(j)==allDays):find(dateBarriers(j+1)==allDays));
            end
            clc
            fprintf('Finding PF for %s to %s\n',datestr(localDays(1)),datestr(localDays(end)))
            
            % let user decide to select date by date or to use all data
            useAllDatesFlag = input(sprintf('Use all data in range(1) or select day by day(0)?: '));
            
            % load data within jth range
            local_u = u(t<=localDays(end)+1);
            local_v = v(t<=localDays(end)+1);
            local_w = w(t<=localDays(end)+1);
            local_direction = direction(t<=localDays(end)+1);
            local_t = t(t<=localDays(end)+1);
            numValues = [];
            
            if useAllDatesFlag
                % find coefficients
                PFinfo.(['cm_',num2str(round(z(ii)*100))]).(['day_',num2str(localDays(1)),'to',num2str(localDays(end))]) =...
                    PF_coefficients([local_u,local_v,local_w,local_direction],bins);
                numValues = numel(local_u);
            else
                numBins = length(bins);
                numDays = length(localDays);
                
                % arrange plots
                figure
                for k = 1:numBins
                    subplot(2,ceil(numBins/2),k)
                    if k ~=numBins
                        title(sprintf('%g < dir < %g',bins(k),bins(k+1)))
                    else
                        title(sprintf('%g < dir < %g',bins(k),bins(1)))
                    end
                    xlabel('Day')
                    ylabel('$^\circ$','interpreter','latex')
                end
                hold all
                changeFig(15)
                
                % initialize cumulative variables
                cumRows = [];
                cumU = [];
                cumV = [];
                cumB = [];
                
                % iterate through days
                for k = 1:numDays
                    
                    % find rows for individual day
                    dayRows = find(localDays(k) == floor(local_t));
                    
                    % find day variables
                    day_u = local_u(dayRows);
                    day_v = local_v(dayRows);
                    day_w = local_w(dayRows);
                    day_direction = local_direction(dayRows);
                    day_t = local_t(dayRows);
                    
                    % find day coefficients
                    in = [day_u, day_v, day_w, day_direction];
                    coef = PF_coefficients(in,bins);
                    
                    % plot day coefficients
                    for m = 1:numBins
                        if m<numBins
                            b1 = coef.(sprintf('degrees_%g_to_%g',bins(m),bins(m+1)))(2);
                            b2 = coef.(sprintf('degrees_%g_to_%g',bins(m),bins(m+1)))(3);
                        else
                            b1 = coef.(sprintf('degrees_%g_to_%g',bins(m),bins(1)))(2);
                            b2 = coef.(sprintf('degrees_%g_to_%g',bins(m),bins(1)))(3);
                        end
                        pitch = asin(-b1/sqrt(1+b1^2))*180/pi;  % put pitch and roll in degrees
                        roll = asin(b2/sqrt(1+b2^2))*180/pi;
                        subplot(2,ceil(numBins/2),m)
                        hold on
                        plot(k,pitch,'x')
                        plot(k,roll,'gd')
                        xlim([0 numDays])
                    end
                    
                    % ask user to approve current day
                    useDayFlag = input('Use current day for cumulative calculation (1=yes, 0=no)?: ');
                    
                    if useDayFlag
                        % find cumulative data
                        cumRows = [cumRows; dayRows];
                        cumU = local_u(cumRows);
                        cumV = local_v(cumRows);
                        cumW = local_w(cumRows);
                        cumDir = local_direction(cumRows);
                        
                        % find cumulative coefficients
                        in = [cumU, cumV, cumW, cumDir];
                        coef = PF_coefficients(in,bins);
                        
                        % plot cumulative coefficients
                        for m = 1:numBins
                            if m<numBins
                                b1Cum = coef.(sprintf('degrees_%g_to_%g',bins(m),bins(m+1)))(2);
                                b2Cum = coef.(sprintf('degrees_%g_to_%g',bins(m),bins(m+1)))(3);
                            else
                                b1Cum = coef.(sprintf('degrees_%g_to_%g',bins(m),bins(1)))(2);
                                b2Cum = coef.(sprintf('degrees_%g_to_%g',bins(m),bins(1)))(3);
                            end
                            pitchCum(k,m) = asin(-b1Cum/sqrt(1+b1Cum^2))*180/pi;  % put pitch and roll in degrees
                            rollCum(k,m) = asin(b2Cum/sqrt(1+b2Cum^2))*180/pi;
                            subplot(2,ceil(numBins/2),m)
                            hold on
                            % find X and Y pitch components.  Delete 0's
                            pitchY = pitchCum(1:k,m);
                            pitchX = 1:k;
                            pitchX(pitchY==0) = [];
                            pitchY(pitchY==0) = [];
                            % plot Pitch
                            h(1) = plot(pitchX,pitchY,'b.-');
                            % find X and Y roll components.  Delete 0's
                            rollY = rollCum(1:k,m);
                            rollX = 1:k;
                            rollX(rollY==0) = [];
                            rollY(rollY==0) = [];
                            % plot roll
                            h(2) = plot(rollX,rollY,'g.-');
                            xlim([0 numDays])
                            legend(h,'Pitch','Roll')
                            changeFig(15,2,8)
                        end
                    else
                        changeFig(15,2,8)
                    end
                end
                numValues = numel(cumU);
                % store PF info
                PFinfo.(['cm_',num2str(round(z(ii)*100))]).(['day_',num2str(localDays(1)),'to',num2str(localDays(end))]) = coef;
            end
        end
    end
    % save PF info
    save([info.rootFolder,filesep,info.siteFolder,filesep,'PFinfo.mat'],'PFinfo')
else % load PF info
    load([info.rootFolder,filesep,info.siteFolder,filesep,'PFinfo.mat'])
end

% display PF info
clc; 
disp('Verify Global Planar Fit Coefficients'); display(' ');
% get heights
heights = fields(PFinfo);

% iterate through heights
for ii = 1:length(heights)
   fprintf('\nSonic Height: %s',heights{ii}); 
   
   % get dates
   dateStrings = fields(PFinfo.(heights{ii}));
   
   % iterate through dates
   for j = 1:length(dateStrings)
       
       % display date range
       fprintf('\t%s to %s\n',datestr(str2double(dateStrings{j}(5:10))),datestr(str2double(dateStrings{j}(13:18))))
       
       % find bins
       bins = fields(PFinfo.(heights{ii}).(dateStrings{j}));
       
       % iterate through bins
       for k = 1:length(bins)
           fprintf('\t\t%s',bins{k})
           b0 = PFinfo.(heights{ii}).(dateStrings{j}).(bins{k})(1);
           b1 = PFinfo.(heights{ii}).(dateStrings{j}).(bins{k})(2);
           b2 = PFinfo.(heights{ii}).(dateStrings{j}).(bins{k})(3);
           pitch = asin(-b1/sqrt(1+b1^2))*180/pi;  % put pitch and roll in degrees
           roll = asin(b2/sqrt(1+b2^2))*180/pi;
           fprintf('\t\t\tb0=%.03g, b1=%.03g, b2=%.03g\n\t\t\tpitch=%.03g deg, roll=%.03g deg\n',b0,b1,b2,pitch,roll)
           PFinfo.infoString{1,ii} = 'Global PF Info';
           PFinfo.infoString{2,ii} = sprintf('Sonic Height: %s m',heights{ii});
           PFinfo.infoString{(j-1)*3+3,ii} = sprintf('%s to %s',datestr(str2double(dateStrings{j}(5:10))),datestr(str2double(dateStrings{j}(13:18))));
           PFinfo.infoString{(j-1)*3+4,ii} = sprintf('%s',bins{k});
           PFinfo.infoString{(j-1)*3+5,ii} = sprintf('b0=%.03g, b1=%.03g, b2=%.03g -- pitch=%.03g deg, roll=%.03g deg',b0,b1,b2,pitch,roll);
       end
   end
end

% ask user to verify information
disp(' ');
PFflag = input('Is this correct (1 for yes, 0 for no): ');

if ~PFflag 
    error('Program stopped by user. Check PF coefficients')
end

clc
flag3 = input('Okay to begin analysis (1 for yes, 0 for no): ');
if ~flag3
    error('Program stopped by user.  Check inputs and restart.')
end
end