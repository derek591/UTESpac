function [rotatedSonicData, output, dataInfo] = sonicRotation(output,data,sensorInfo,info,dataInfo,tableNames,PFinfo)
% sonicRotation applies a single sector, local (single .csv file) or global, multi-sector (all .csv files in siteFolder)
% planar fit and yaw rotation to the sonic data.  The global coefficients are created in PFinfo.m (line ?? in UTESpac)
% and the global option is controlled by info.PF.recalculateGlobalCoefficients (line ?? in UTESpac). 
if isfield(sensorInfo,'u')
    fprintf('\nPlanar Fitting and Performing Yaw Rotation on Sonics\n')
    % find number of sonics
    numSonics = size(sensorInfo.u,1);
        
    % iterate through all sonics
    for ii = 1:numSonics
        try
            % find sonic information
            tble = sensorInfo.u(ii,1);
            sonHeight = sensorInfo.u(ii,3);
            uCol = sensorInfo.u(sensorInfo.u(:,3)==sonHeight,2);
            vCol = sensorInfo.v(sensorInfo.v(:,3)==sonHeight,2);
            wCol = sensorInfo.w(sensorInfo.v(:,3)==sonHeight,2);
            
            % find unaveraged sonic values to be rotated
            u = data{1,tble}(:,uCol);
            v = data{1,tble}(:,vCol);
            w = data{1,tble}(:,wCol);
            t = data{1,tble}(:,1);
            
            % find averaged direction and time
            dirAvg = output.spdAndDir(:,strcmp(output.spdAndDirHeader,[num2str(sonHeight),'m direction']));
            tAvg = output.spdAndDir(:,1);
            
            % create PF wind to store rotated wind date from ith sonic
            windPF = nan(size(u,1),3);
            %------------------------------------ PLANAR FIT (wbar ->0) -----------------------------------------------------------
            
            % --- GLOBAL PF CALCULATION
            if strcmp(info.PF.globalCalculation,'global') && exist('PFinfo','var') && isfield(PFinfo,sprintf('cm_%g',sonHeight*100))
                
                % put global pf info into data info column
                if ii == 1
                    dataInfo(1:size(PFinfo.infoString,1),end+1:end+size(PFinfo.infoString,2)) = PFinfo.infoString;
                end
                
                % find all date bins
                dateBins = fields(PFinfo.(sprintf('cm_%g',sonHeight*100)));
                
                % iterate through datebins
                for jj = 1:length(dateBins)
                    
                    % extract date numbers from field name.  Replace '_' and 'to' from bin name to allow sscanf
                    jth_PF_date_bin = sscanf(strrep(strrep(dateBins{jj},'_',' '),'to',' '),'%*s %d %d');
                    local_dates_span = [min(tAvg), max(tAvg)];
                    
                    % check to see if local data falls within date bin.  Otherwise, continue.
                    if local_dates_span(1)>=jth_PF_date_bin(1) && local_dates_span(1)<=jth_PF_date_bin(2) || local_dates_span(2)>=jth_PF_date_bin(1) && local_dates_span(2)<=jth_PF_date_bin(2)
                        
                        % find all direcition bins
                        dirBins = fields(PFinfo.(sprintf('cm_%g',sonHeight*100)).(dateBins{jj}));
                        
                        % iterate though all direction bins
                        for kk = 1:length(dirBins)
                            
                            % find local direction bin
                            kth_dir_bin = sscanf(strrep(dirBins{kk},'_',' '),'%*s %d %*s %d');
                            
                            % find local coefficients
                            localCoef = PFinfo.(sprintf('cm_%g',sonHeight*100)).(dateBins{jj}).(dirBins{kk});
                            b1 = localCoef(1);
                            b2 = localCoef(2);
                            
                            % find Values of B C and D matrices
                            p31 = -b1/(b1^2+b2^2+1)^0.5;
                            p32 = -b2/(b1^2+b2^2+1)^0.5;
                            p33 = 1/(b1^2+b2^2+1)^0.5;
                            
                            % Wilczak Large Angle
                            sin_alpha = p31;
                            cos_alpha = (p32^2+p33^2)^0.5;
                            sin_beta = -p32/(p32^2+p33^2)^0.5;
                            cos_beta = p33/(p32^2+p33^2)^0.5;
                            
                            % populate C and D matrices
                            C = [1 0 0; 0 cos_beta -sin_beta; 0 sin_beta cos_beta];
                            D = [cos_alpha 0 sin_alpha; 0 1 0; -sin_alpha 0 cos_alpha];
                            
                            % find P
                            P = D'*C';
                                                       
                            % expand mean wind directions to turbulence data
                            % (http://www.mathworks.com/matlabcentral/newsreader/view_thread/44023)
                            repititions = length(t)/length(tAvg)*ones(length(tAvg),1);
                            h = [];
                            h(cumsum(repititions))=1;
                            mean_dir_20Hz=dirAvg(cumsum(h)-h+1,:);
                            
                            % find rows within jjth date bin
                            dateLogical = zeros(size(t));
                            dateLogical(t>=jth_PF_date_bin(1) & t<=jth_PF_date_bin(2)) = 1;
                            
                            % find rows within kkth direction bin
                            dirLogical = zeros(size(t));
                            if kth_dir_bin(1) - kth_dir_bin(2) == 0  % one sector
                                dirLogical = 1;
                            elseif kth_dir_bin(2) > kth_dir_bin(1) % bin does not include North
                                dirLogical(mean_dir_20Hz>kth_dir_bin(1) & mean_dir_20Hz<kth_dir_bin(2)) = 1;
                            else % bin includes North
                                dirLogical(mean_dir_20Hz>kth_dir_bin(1) | mean_dir_20Hz<kth_dir_bin(2)) = 1;
                            end
                            
                            % combine dateLogical and dirLogical
                            rowLogical = dateLogical+dirLogical;
                            rowLogical(rowLogical<2) = 0;
                            rowLogical = logical(rowLogical);
                            
                            % apply local planar fit to data
                            windPF(rowLogical,:) = (P*[u(rowLogical) v(rowLogical) w(rowLogical)]')';
                        end
                    else % if no local data falls within jjth date bin, continue
                        continue
                    end
                end
                % failed global calculation due to non-existant PFinfo
            elseif strcmp(info.PF.globalCalculation,'global') && info.PF.recalculateGlobalCoefficients && ~exist('PFinfo','var')
                error('PFinfo variable failed to load.  Check/re-run findGlobalPR.mat')
                
                % failed global calculation due to non-existant PFinfo at ith height
            elseif strcmp(info.PF.globalCalculation,'global') && info.PF.recalculateGlobalCoefficients && exist('PFinfo','var') && ~isfield(PFinfo,sprintf('cm_%g',sonHeight*100))
                error('PF info at %g m does not exist.  Check/re-run findGlobalPR.mat',sonHeight)
                
                % --- LOCAL PF CALCULATION
            else
                % initialize dataInfo column to store local coefficients in
                if ii==1 || ~exist('dataInfoCol','var')
                    dataInfoCol = size(dataInfo,2)+1;
                    dataInfo{1,dataInfoCol} = 'Local PF Coefficients';
                end
                
                % find spike and nan flags.  combine flags for u, v and w with logical()
                spikeFlag = logical(output.([tableNames{tble},'SpikeFlag'])(:,uCol) + ...
                    output.([tableNames{tble},'SpikeFlag'])(:,vCol) + ...
                    output.([tableNames{tble},'SpikeFlag'])(:,wCol));
                nanFlag = logical(output.([tableNames{tble},'NanFlag'])(:,uCol) + ...
                    output.([tableNames{tble},'NanFlag'])(:,vCol) + ...
                    output.([tableNames{tble},'NanFlag'])(:,wCol));
                
                % find wind flag
                windFlagName = sprintf('%gm flag',sonHeight); % name of wind flag column
                windFlagColumn = strncmp(output.spdAndDirHeader,windFlagName,length(windFlagName)); % wind flag column
                windFlag = output.spdAndDir(:,windFlagColumn);
                
                % combine all flags
                FLAG = logical(spikeFlag+nanFlag+windFlag);
                
                % find averaged sonic values to compute planar fit.  Columns = ~FLAG
                ubar = output.(tableNames{tble})(~FLAG,uCol);
                vbar = output.(tableNames{tble})(~FLAG,vCol);
                wbar = output.(tableNames{tble})(~FLAG,wCol);
                
                % find number of records
                flen = length(ubar);
                
                % find values to populate matrix equation
                su = sum(ubar);
                sv = sum(vbar);
                sw = sum(wbar);
                suv = sum(ubar.*vbar);
                suw = sum(ubar.*wbar);
                svw = sum(vbar.*wbar);
                su2 = sum(ubar.*ubar);
                sv2 = sum(vbar.*vbar);
                
                % [H](b) = (g)
                H = [flen su sv; su su2 suv; sv suv sv2];
                g = [sw suw svw]';
                
                % solve for b coefficients 
                coef = linsolve(H,g);
                
                % only b1 and b2 are necessary for planar fit
                b1 = coef(2);
                b2 = coef(3);
                
                % output pitch and roll angles to command window and dataInfo
                pitch = asin(-b1/sqrt(1+b1^2))*180/pi;  
                roll = asin(b2/sqrt(1+b2^2))*180/pi;
                fprintf('Sonic @ %gm\t pitch=%.3g\t roll=%.3g degrees.\n',sensorInfo.u(ii,3),pitch,roll);
                dataInfo{ii+1,dataInfoCol} = sprintf('%gm\t b0=%.3g\t b1=%.3g\t b2=%.3g -\t pitch=%.3g\t roll=%.3g deg',sonHeight,coef(1),coef(2),coef(3),pitch,roll);
                              
                % find Values of B C and D matrices
                p31 = -b1/(b1^2+b2^2+1)^0.5;
                p32 = -b2/(b1^2+b2^2+1)^0.5;
                p33 = 1/(b1^2+b2^2+1)^0.5;
                
                % Wilczak Large Angle
                sin_alpha = p31;
                cos_alpha = (p32^2+p33^2)^0.5;
                sin_beta = -p32/(p32^2+p33^2)^0.5;
                cos_beta = p33/(p32^2+p33^2)^0.5;
                
                % populate C and D matrices
                C = [1 0 0; 0 cos_beta -sin_beta; 0 sin_beta cos_beta];
                D = [cos_alpha 0 sin_alpha; 0 1 0; -sin_alpha 0 cos_alpha];
                
                % find P
                P = D'*C';
                
                % apply local planar fit to data
                windPF = (P*[u v w]')';
            end
            %------------------------------------ YAW ROTATION (vbar ->0) ----------------------------------------------------------
            
            % find averaged PF data
            windPFavg = simpleAvg([windPF(:,1:2) t],info.avgPer);
            
            % partition averaged PF sonic data
            uPFbar = windPFavg(:,1);
            vPFbar = windPFavg(:,2);
            
            % find cos and sin of gamma
            cosGamma = uPFbar./(uPFbar.^2+vPFbar.^2).^0.5;
            sinGamma = vPFbar./(uPFbar.^2+vPFbar.^2).^0.5;
            
            % define time step
            dt = floor(t(1,end)) + datenum(0,0,0,0,info.avgPer,0);
            
            % initialize lower index
            I0 = 1;
            
            % find number of averaging periods
            N = (ceil(t(end))-floor(t(1)))/(info.avgPer/(24*60));
            
            % initialize rotated data
            if ii == 1
                rotatedSonicData = nan(size(windPF));
            else
                rotatedSonicData(:,3*ii-2:3*ii) = nan(size(windPF,1),3);
            end
            % create output header
            output.rotatedSonicHeader{1,3*ii-2} = sprintf('%gm:u',sonHeight);
            output.rotatedSonicHeader{1,3*ii-1} = sprintf('%gm:v',sonHeight);
            output.rotatedSonicHeader{1,3*ii} = sprintf('%gm:w',sonHeight);
            
            % iterate through all rows
            for row = 1:N
                % find index of last value in the averaging period
                I1 = find(t >= dt,1,'first') - 1;
                
                % check to make sure an I1 is found.  This is often not the case over the last averaging period
                if isempty(I1) == true
                    I1 = size(windPF,1);
                end
                
                % find rotation matrix
                M = [cosGamma(row) sinGamma(row) 0; -sinGamma(row) cosGamma(row) 0; 0 0 1];
                
                % rotate data over averaging period
                rotatedSonicData(I0:I1,3*ii-2:3*ii) = (M*windPF(I0:I1,:)')';
                
                % set up next iteration
                I0 = I1+1;
                dt = dt+datenum(0,0,0,0,info.avgPer,0);
            end
        catch err
            message = strcat(err.message,'@ line ',num2str(err.stack.line),' UNABLE TO APPLY SONIC ROTATIONS AT ',num2str(sonHeight),'m');
            warning(message)
            if isempty(output.warnings{1})
                output.warnings{1,1} = message;
            else
                output.warnings{end+1,1} = message;
            end
        end
    end
    % average and store rotated sonic data
    try
        temp = simpleAvg([rotatedSonicData t],info.avgPer);
        output.rotatedSonic = temp(:,1:end-1);
    catch err
        message = strcat(err.message,'@ line ',num2str(err.stack.line),' UNABLE TO APPLY SONIC ROTATIONS AT All HEIGHTS');
        warning(message)
        rotatedSonicData = [];
        if isempty(output.warnings{1})
            output.warnings{1,1} = message;
        else output.warnings{end+1,1} = message;
        end
        pause(3);
    end
else
    rotatedSonicData = [];
end
end