function [output, raw] = fluxes(data, rotatedSonicData, info, output, sensorInfo,tableNames)
% fluxes computes all turbulent statistics.  Spike and NaN flags are used from sonic and finewires to nan-out flagged
% statistics.
fprintf('\nComputing Fluxes\n')
% ensure that sonics exist
if ~isfield(sensorInfo,'u')
    raw = [];
    return
end
try
    %---------------- FIND REFERENCE VALUES
    % constants
    Rd = 287.058;  % [J/K/kg] Gas constant for air
    Rv = 461.495;  % [J/K/kg] Gas constant for water vapor
    
    % find sonic time stamps
    t = data{1,sensorInfo.u(1,1)}(:,1); % time stamps for sonics
    
    % find reference pressure
    zRef = min(sensorInfo.u(:,3));  % zRef is lowest sonic level
    altitude = info.siteElevation;
    Pref = 101325*(1-2.25577*10^-5*(altitude+zRef))^5.25588/1000; % find pRef from elevation: http://www.engineeringtoolbox.com/air-altitude-pressure-d_462.html
    
    % find pressure, use calculated Pref if barometer does not exist, is all nans or deviates from Pref by more than 5%
    if isfield(sensorInfo,'P')
        zRef = sensorInfo.P(1,3);   % Pressure Height
        Ptble = sensorInfo.P(1,1);  % Pressure table
        Pcol = sensorInfo.P(1,2);   % Pressure Column
        PkPa = data{1,Ptble}(:,Pcol);  % Pressure in kPa or mbar
        Ptime = data{1,Ptble}(:,1); % pressure time stamps
        
        % ensure that Pressure is in kPa
        if nanmedian(PkPa)>200
            PkPa = PkPa./10;
        end
        
        % store PkPa in raw tables
        if info.saveRawConditionedData
            raw.P = [Ptime, PkPa]; % u
        end
        
        % average P down to avgPer
        PkPaAvg = simpleAvg([PkPa, Ptime],info.avgPer,0);
        clear P Ptime
        if sum(isnan(PkPaAvg))/numel(PkPaAvg) < 1
            fprintf('Barometer found. Median Pressure = %0.03g kPa\n',nanmedian(PkPa))
        end
    end
    
    if ~isfield(sensorInfo,'P') || sum(isnan(PkPaAvg))/numel(PkPaAvg) == 1 || abs((nanmedian(PkPaAvg) - Pref)/Pref) > 0.05
        zRef = min(sensorInfo.u(:,3));  % zRef is lowest sonic level
        altitude = info.siteElevation;
        Pref = 101325*(1-2.25577*10^-5*(altitude+zRef))^5.25588/1000; % find pRef from elevation: http://www.engineeringtoolbox.com/air-altitude-pressure-d_462.html
        PkPaAvg = Pref*ones(size(output.rotatedSonic(:,1))); % create pressure array
        fprintf('No barometer found or barometer values all NaNs or median barometer values deviate from reference pressure by more than 5%%.\nReference Pressure of %0.03g kPa calculated from site elevation (defined in siteInfo.m) will be used.\n',Pref)
    end
    
    % find reference temperature nearest to zRef
    if isfield(sensorInfo,'T')  % look for slow repsonse temperature first (HMP45/155)
        % find temperature closest to zRef
        [~, Tindex] = min(abs(sensorInfo.T(:,3)-zRef));  % find T closest to zRef
        T_tble = sensorInfo.T(Tindex,1);  % Find sonic corresponding to pressure level
        T_col = sensorInfo.T(Tindex,2);
        Tref_K = data{1,T_tble}(:,T_col); % reference sonic temperature in C or K
        Tref_time = data{1,T_tble}(:,1);
        Tref_Kavg = simpleAvg([Tref_K, Tref_time],info.avgPer,0); % average down to avgPer
        
        
        % ensure Tref is in K
        if nanmedian(Tref_Kavg)<200
            Tref_Kavg = Tref_Kavg + 273.15;  % put sonTref_K in K
        end
        
        if sum(isnan(Tref_Kavg))/numel(Tref_Kavg) < 1
            fprintf('Slow Repsonse Temperature found. Median Reference Temperature = %0.03g C\n',nanmedian(Tref_Kavg)-273.15)
        end
    end
    
    if ~isfield(sensorInfo,'T') || sum(isnan(Tref_Kavg))/numel(Tref_Kavg) == 1 % if slow response T unavailable, uses closest sonic to zRef
        % find sonic reference temperature closest to zRef
        [~, Tindex] = min(abs(sensorInfo.u(:,3)-zRef));  % find sonic closest to zRef
        T_tble = sensorInfo.Tson(Tindex,1);
        T_col = sensorInfo.Tson(Tindex,2);
        Tref_K = data{1,T_tble}(:,T_col); % reference sonic temperature in C or K
        
        Tref_Kavg = simpleAvg([Tref_K, t],info.avgPer,0); % average down to avgPer
        
        
        % ensure Tref is in K
        if nanmedian(Tref_Kavg)<200
            Tref_Kavg = Tref_Kavg + 273.15;  % put sonTref_K in K
        end
        fprintf('No Slow Repsonse Temperature found, or slow response temperature is all NaNs! Median Reference Temperature from Sonic = %0.03g C\n',nanmedian(Tref_Kavg)-273.15)
    end
    
    % find reference specific humidity, this will be updated at each level with the nearest HMP if available
    if isfield(sensorInfo,'RH')  % try slow response RH sensor first
        
        % find closes RH measurement to zRef
        [~, RHindex] = min(abs(sensorInfo.RH(:,3)-zRef));
        RHtble = sensorInfo.RH(RHindex,1);
        RHcol = sensorInfo.RH(RHindex,2);
        RH = data{1,RHtble}(:,RHcol);  % RH closest to zRef
        RHtime = data{1,RHtble}(:,1);
        
        % interpolate RH to same time stamps as sonic_time
        RHavg = simpleAvg([RH, RHtime],info.avgPer,1);
        RHavg_t = RHavg(:,2);
        RHavg(:,2) = [];
        clear RH
        
        % find average specific humidity
        qRefavg = RHtoSpecHum(RHavg,PkPaAvg,Tref_Kavg); % [kg/kg]
        
        % ensure that qRefavg is not all NaNs
        if sum(isnan(qRefavg))<numel(qRefavg)
            % check output here: http://www.rotronic.com/humidity_measurement-feuchtemessung-mesure_de_l_humidite/humidity-calculator-feuchterechner-mr
            % interpolate to sonic time stamps.  Pad qRefavg at beginning to allow interpolation at all time stamps
            x = RHavg_t(~isnan(qRefavg));
            x = [floor(x(1)); x];
            y = qRefavg(~isnan(qRefavg));
            y = [y(1); y];
            qRefFast = interp1(x,y,t);  % interpolate qRef to Sonic Frequency
            fprintf('Slow Repsonse RH found! Median Reference Specific Humidity = %0.03g g/kg\n',1000*nanmedian(qRefavg))
        end
    else
        qRefavg = nan;
    end
    
    % find qref from IRGA or KH2O
    if ~isfield(sensorInfo,'RH') && (isfield(sensorInfo,'irgaH2O') || isfield(sensorInfo,'KH2O')) 
        
        if isfield(sensorInfo,'irgaH2O')
            rhoH2O_tble = sensorInfo.irgaH2O(1,1);
            rhoH2O_col = sensorInfo.irgaH2O(1,2);
            rhovIRGA = data{1,rhoH2O_tble}(:,rhoH2O_col);
            rhov_t = data{1,rhoH2O_tble}(:,1);
            rhovIRGAavg = simpleAvg([rhovIRGA,rhov_t],info.avgPer);
            disp('No slow-response RH found. qRef calculated from EC150')
        elseif isfield(sensorInfo,'KH2O')
            rhoH2O_tble = sensorInfo.KH2O(1,1);
            rhoH2O_col = sensorInfo.KH2O(1,2);
            rhovIRGA = data{1,rhoH2O_tble}(:,rhoH2O_col);
            rhov_t = data{1,rhoH2O_tble}(:,1);
            rhovIRGAavg = simpleAvg([rhovIRGA,rhov_t],info.avgPer);
            disp('No slow-response RH found. qRef calculated from KH2O')
        elseif isfield(sensorInfo,'LiH2O')
            rhoH2O_tble = sensorInfo.LiH2O(1,1);
            rhoH2O_col = sensorInfo.LiH2O(1,2);
            rhovIRGA = data{1,rhoH2O_tble}(:,rhoH2O_col)*18/1000;  % Multiply by 18/1000 to go from mmol/mol to g/m^3
            rhov_t = data{1,rhoH2O_tble}(:,1);
            rhovIRGAavg = simpleAvg([rhovIRGA,rhov_t],info.avgPer);
            disp('No slow-response RH found. qRef calculated from EC150') 
        end
        qRefavg = (rhovIRGAavg(:,1)./1000)./(PkPaAvg*1000./(Rd*Tref_Kavg)); % kg/kg  rhoH2O/rhoAir
        x = rhovIRGAavg(~isnan(qRefavg),2);
        x = [floor(x(1)); x];
        y = qRefavg(~isnan(qRefavg));
        y = [y(1); y];
        qRefFast = interp1(x,y,t);  % interpolate qRef to Sonic Frequency
        fprintf('Median Reference Specific Humidity = %0.03g g/kg\n',1000*nanmedian(qRefavg))
    elseif ~isfield(sensorInfo,'RH') ||  sum(isnan(qRefavg))/numel(qRefavg)==1 % if no humidity measurements exist, use info.qRef
        qRefavg = info.qRef./1000*ones(size(Tref_Kavg));
        qRefFast = info.qRef./1000*ones(size(t));
        fprintf('No slow-response nor fast response RH found (or humidity measurement is all NaNs). qRef = %0.02g g/kg defined in INORMATION block is being used.',info.qRef)
    end
    % store specific humidity in output
    
    % find moist-air density, dry-air density and vapor density
    PvAvg = qRefavg.*PkPaAvg./0.622;  % partial pressure of water vapor (kPa)
    PdAvg = PkPaAvg-PvAvg; % partial pressure of dry air
    rhodAvg = 1000*PdAvg./(Rd*Tref_Kavg); % density of dry air (kg/m^3)
    rhovAvg = 1000*PvAvg./(Rv*Tref_Kavg); % density of water vapor (kg/m^3)
    rhoAvg = rhodAvg+rhovAvg; % density of moist air (kg/m^3)
    fprintf('Calculated Moist Air Density = %0.03g kg/m^3, Dry Air Density = %0.03g kg/m^3\n',nanmedian(rhoAvg),nanmedian(rhodAvg))
    
    
    %------------- ITERATE THROUGH SONICS
    % find number of sonics
    numSonics = size(sensorInfo.u,1);
    
    for ii = 1:numSonics
        if ii==1
            
            % find total number of averaging periods
            N = round((t(end)-t(1))/(info.avgPer/(24*60)));
            
            % find breakpoints for detrending.  size(bp,1) = size(N,1) + 1
            bp = round(linspace(0,numel(t),N+1));
            
            % initialize flux matrices
            H = nan(N,1);  % kinematic sensible heat flux
            Hlat = nan(N,1); % kinematic sensible, lateral heat flux
            tau = nan(N,1); % momentum flux
            tke = nan(N,1); % turbulent kinetic energy
            LHflux = nan(N,1); % latent heat flux
            CO2flux = nan(N,1); % CO2 flux
            derivedT = nan(N,1);  % matrix for derived temperatures
            L = nan(N,1); % Obukhov Length
            sigma = nan(N,1); % standard deviations (u, v, w)
            
            % initialize raw variable matrices
            if info.saveRawConditionedData
                raw.u = nan(bp(end),numSonics); % u
                raw.v = nan(bp(end),numSonics); % v
                raw.w = nan(bp(end),numSonics); % w
                raw.uPF = nan(bp(end),numSonics);  % u planar fit and yaw corrected
                raw.vPF = nan(bp(end),numSonics);  % v planar fit and yaw corrected
                raw.wPF = nan(bp(end),numSonics);  % w planar fit and yaw corrected
                raw.sonTs = nan(bp(end),numSonics); % Temperature from Sonic
                raw.uPF_Prime = nan(bp(end),numSonics);  % u' planar fit and yaw corrected
                raw.vPF_Prime = nan(bp(end),numSonics); % v' planar fit and yaw corrected
                raw.wPF_Prime = nan(bp(end),numSonics); % w' planar fit and yaw corrected
                raw.sonTsPrime = nan(bp(end),numSonics); % Ts' from sonic
                raw.t = t; % serial time stamp
                raw.z = nan(1,numSonics); % sonic heights
                
                % finewires
                if isfield(sensorInfo,'fw')
                    raw.fwTh = nan(bp(end),numSonics); % theta from FW
                    raw.fwT = nan(bp(end),numSonics); % T from FW
                    faw.fwTPrime = nan(bp(end),numSonics); % T' from FW
                    raw.fwThPrime = nan(bp(end),numSonics); % theta' from finewire
                end
                
                
                % water vapor
                if isfield(sensorInfo,'irgaH2O')
                    numH2OSensors = size(sensorInfo.irgaH2O,1);
                    raw.rhov = nan(bp(end),numH2OSensors); % H2O
                    raw.rhovPrime = nan(bp(end),numH2OSensors); % H2O'
                elseif isfield(sensorInfo,'LiH2O')
                    numH2OSensors = size(sensorInfo.LiH2O,1);
                    raw.rhov = nan(bp(end),numH2OSensors); % H2O
                    raw.rhovPrime = nan(bp(end),numH2OSensors); % H2O'
                elseif isfield(sensorInfo,'KH2O')
                    numH2OSensors = size(sensorInfo.KH2O,1);
                    raw.rhov = nan(bp(end),numH2OSensors); % H2O
                    raw.rhovPrime = nan(bp(end),numH2OSensors); % H2O'
                end
                
                % CO2
                if isfield(sensorInfo,'irgaCO2')
                    numCO2Sensors = size(sensorInfo.irgaCO2,1);
                    raw.rhoCO2 = nan(bp(end),numCO2Sensors); % CO2
                    raw.rhoCO2Prime = nan(bp(end),numCO2Sensors); % CO2'
                elseif isfield(sensorInfo,'LiCO2')
                    numCO2Sensors = size(sensorInfo.LiCO2,1);
                    raw.rhoCO2 = nan(bp(end),numCO2Sensors); % CO2
                    raw.rhoCO2Prime = nan(bp(end),numCO2Sensors); % CO2'
                end
                
            else
                raw = [];
            end
            
            % initialize headers
            Hheader = cell(1);
            HlatHeader = cell(1);
            tauHeader = cell(1);
            tkeHeader = cell(1);
            LHfluxHeader = cell(1);
            CO2fluxHeader = cell(1);
            derivedTheader = cell(1);
            Lheader = cell(1);
            sigmaHeader = cell(1);
        end
        try
            
            % find sonic information
            tble = sensorInfo.u(ii,1);
            sonHeight = sensorInfo.u(ii,3);
            uCol = sensorInfo.u(sensorInfo.u(:,3)==sonHeight,2);
            vCol = sensorInfo.v(sensorInfo.v(:,3)==sonHeight,2);
            wCol = sensorInfo.w(sensorInfo.v(:,3)==sonHeight,2);
            TsCol = sensorInfo.Tson(sensorInfo.Tson(:,3)==sonHeight,2);
            
            % find sonic flag information (averaged values!)
            nanFlagTableName = [tableNames{tble},'NanFlag'];
            spikeFlagTableName = [tableNames{tble},'SpikeFlag'];
            uNanFlag = output.(nanFlagTableName)(:,uCol);
            uSpikeFlag = output.(spikeFlagTableName)(:,uCol);
            vNanFlag = output.(nanFlagTableName)(:,vCol);
            vSpikeFlag = output.(spikeFlagTableName)(:,vCol);
            wNanFlag = output.(nanFlagTableName)(:,wCol);
            wSpikeFlag = output.(spikeFlagTableName)(:,wCol);
            TsNanFlag = output.(nanFlagTableName)(:,TsCol);
            TsSpikeFlag = output.(nanFlagTableName)(:,TsCol);
            
            % check for diagnostic
            if isfield(sensorInfo,'sonDiagnostic')
                sonicDiagnosticCol = sensorInfo.sonDiagnostic(sensorInfo.sonDiagnostic(:,3)==sonHeight,2);
                if isempty(sonicDiagnosticCol)
                    sonicDiagnosticFlag = zeros(size(uNanFlag));
                else
                    sonicDiagnosticFlag = output.(tableNames{tble})(:,sonicDiagnosticCol);
                    sonicDiagnosticFlag(sonicDiagnosticFlag < info.diagnosticTest.meanSonicDiagnosticLimit) = 0;
                end
            else
                sonicDiagnosticFlag = zeros(size(uNanFlag));
            end
            sonicDiagnosticFlag(isnan(sonicDiagnosticFlag)) = 0;
            
            % find unrotated sonic values
            u = data{1,tble}(:,uCol);
            v = data{1,tble}(:,vCol);
            w = data{1,tble}(:,wCol);
            unrotatedSonFlag = logical(wNanFlag+wSpikeFlag+sonicDiagnosticFlag); % total sonic flag for unrotated calculations
            
            % find rotated sonic columns
            uCol = strcmp(output.rotatedSonicHeader,strcat(num2str(sonHeight),'m:u'));
            vCol = strcmp(output.rotatedSonicHeader,strcat(num2str(sonHeight),'m:v'));
            wCol = strcmp(output.rotatedSonicHeader,strcat(num2str(sonHeight),'m:w'));
            
            % load rotated sonic values
            uPF = rotatedSonicData(:,uCol);
            vPF = rotatedSonicData(:,vCol);
            wPF = rotatedSonicData(:,wCol);
            rotatedSonFlag = logical(uNanFlag+uSpikeFlag+vNanFlag+vSpikeFlag+wNanFlag+wSpikeFlag+sonicDiagnosticFlag); % total sonic flag for rotated calculations
            
            % find sonic temperature
            Tson = data{1,tble}(:,TsCol);
            if nanmedian(Tson) > 250  % ensure temperature is in C
                Tson = Tson - 273.15;
            end
            TsonFlag = logical(TsNanFlag+TsSpikeFlag);
            
            % find fw
            if isfield(sensorInfo,'fw')
                fwCol = sensorInfo.fw(sensorInfo.fw(:,3)==sonHeight,2);
                fw = data{1,tble}(:,fwCol);
                nanFlagTableName = [tableNames{tble},'NanFlag'];
                spikeFlagTableName = [tableNames{tble},'SpikeFlag'];
                fwNanFlag = output.(nanFlagTableName)(:,fwCol);
                fwSpikeFlag = output.(spikeFlagTableName)(:,fwCol);
                fwFlag = logical(fwNanFlag+fwSpikeFlag);
            else
                fw = [];
                fwFlag = [];
            end
            
            % find specific humidity at level if HMP exists
            if isfield(sensorInfo,'RH') && isfield(sensorInfo,'T') && any(sensorInfo.RH(:,3) == sonHeight)
                
                % find RH at sensor level
                RHtble = sensorInfo.RH(sensorInfo.RH(:,3)==sonHeight,1);
                RHcol = sensorInfo.RH(sensorInfo.RH(:,3)==sonHeight,2);
                RHavg = simpleAvg([data{1,RHtble}(:,RHcol), data{1,RHtble}(:,1)],info.avgPer,0);
                
                % find T at sensor level
                Ttble = sensorInfo.T(sensorInfo.T(:,3)==sonHeight,1);
                Tcol = sensorInfo.T(sensorInfo.T(:,3)==sonHeight,2);
                Tavg = simpleAvg([data{1,Ttble}(:,Tcol), data{1,Ttble}(:,1)],info.avgPer,0);
                
                % put Tavg in K
                if nanmedian(Tavg) < 250; Tavg = Tavg + 273.15; end
                
                % intialize specific humidity output, populate col 1 with time stamps
                if ~isfield(output,'specificHum')
                    output.specificHum = t(bp(2:end));
                    output.specificHumHeader = cell(1);
                    output.specificHumHeader{1} = 'time';
                end
                
                % find qRefLocal and qRefFastLocal
                qRefavgLocal = RHtoSpecHum(RHavg,PkPaAvg,Tavg);
                x = RHavg_t(~isnan(qRefavgLocal));
                if isempty(x)  % if all NaNs, use qRef (non-local)
                    qRefFastLocal = qRefFast;
                    
                    % store nan'd data at height
                    output.specificHum(:,end+1) = qRefavgLocal*nan;
                    output.specificHumHeader{end+1} = sprintf('%g m: q(g/g)',sonHeight);
                else
                    x = [floor(x(1)); x];
                    y = qRefavgLocal(~isnan(qRefavgLocal));
                    y = [y(1); y];
                    qRefFastLocal = interp1(x,y,t);  % interpolate qRef to Sonic Frequency
                    
                    % store data
                    output.specificHum(:,end+1) = qRefavgLocal;
                    output.specificHumHeader{end+1} = sprintf('%g m: q(g/g)',sonHeight);
                end
                
                % use qRefFastLocal if not all NaNs to find virt temp
                if sum(isnan(qRefFastLocal))/numel(qRefFastLocal) == 1
                    qRefFastLocal = qRefFast;
                end
            else
                qRefFastLocal = qRefFast;
            end
            
            % find fw pot temp
            Gamma = 0.0098; % Dry Lapse Rate. K/m
            if ~isempty(fw)
                thetaFw = fw + Gamma*(sonHeight - zRef);
                temp = simpleAvg([thetaFw t],info.avgPer);
                derivedT(:,end+1) = temp(:,1);
                derivedTheader{end+1} = strcat(num2str(sonHeight),' m: theta_fw');
            else
                thetaFw = [];
            end
            
            % find sonic pot temp
            thetaSon = Tson + Gamma*(sonHeight - zRef);
            temp = simpleAvg([thetaSon t],info.avgPer);
            derivedT(:,end+1) = temp(:,1);
            derivedTheader{end+1} = strcat(num2str(sonHeight),' m: theta_s_son');
            
            % find fw virt, pot temp with qRef
            if ~isempty(thetaFw)
                VthetaFw = thetaFw.*(1+0.61*qRefFastLocal); % stull pg 7
                temp = simpleAvg([VthetaFw t],info.avgPer);
                derivedT(:,end+1) = temp(:,1);
                derivedTheader{end+1} = strcat(num2str(sonHeight),' m: theta_v_fw');
            else
                VthetaFw = [];
            end
            
            % find pot temp from sonic
            thetaSonAir = thetaSon./(1+0.51*qRefFastLocal); % stull pg 7
            temp = simpleAvg([thetaSonAir t],info.avgPer);
            derivedT(:,end+1) = temp(:,1);
            derivedTheader{end+1} = strcat(num2str(sonHeight),' m: theta_son');
            
            % find H2O and CO2 columns if they exist
            if isfield(sensorInfo,'irgaH2O') && ~isempty(sensorInfo.irgaH2O(sensorInfo.irgaH2O(:,3)==sonHeight,2)) % EC150
                irgaH2Ocol = sensorInfo.irgaH2O(sensorInfo.irgaH2O(:,3)==sonHeight,2);
                irgaGasDiagCol = sensorInfo.irgaGasDiag(sensorInfo.irgaGasDiag(:,3)==sonHeight,2);
                irgaH2OSigCol = sensorInfo.irgaH2OsigStrength(sensorInfo.irgaH2OsigStrength(:,3)==sonHeight,2);
                irgaH2O = data{1,tble}(:,irgaH2Ocol);
                nanFlagTableName = [tableNames{tble},'NanFlag'];
                spikeFlagTableName = [tableNames{tble},'SpikeFlag'];
                H2ONanFlag = output.(nanFlagTableName)(:,irgaH2Ocol);
                H2OSpikeFlag = output.(spikeFlagTableName)(:,irgaH2Ocol);
                H2OsigFlag = output.(tableNames{tble})(:,irgaH2OSigCol);
                H2OsigFlag(H2OsigFlag > info.diagnosticTest.H2OminSignal) = 0;
                H2OsigFlag(isnan(H2OsigFlag)) = 0;
                gasDiagFlag = output.(tableNames{tble})(:,irgaGasDiagCol);
                gasDiagFlag(gasDiagFlag < info.diagnosticTest.meanGasDiagnosticLimit) = 0;
                gasDiagFlag(isnan(gasDiagFlag)) = 0;
                H2OFlag = logical(H2ONanFlag+H2OSpikeFlag+H2OsigFlag+gasDiagFlag);
            elseif isfield(sensorInfo,'LiH2O') && ~isempty(sensorInfo.LiH2O(sensorInfo.LiH2O(:,3)==sonHeight,2)) % Li7500
                irgaH2Ocol = sensorInfo.LiH2O(sensorInfo.LiH2O(:,3)==sonHeight,2);
                irgaH2O = data{1,tble}(:,irgaH2Ocol)*0.018; % convert from mmol/m^3 to g/m^3
                nanFlagTableName = [tableNames{tble},'NanFlag'];
                spikeFlagTableName = [tableNames{tble},'SpikeFlag'];
                H2ONanFlag = output.(nanFlagTableName)(:,irgaH2Ocol);
                H2OSpikeFlag = output.(spikeFlagTableName)(:,irgaH2Ocol);
                % check for LiCor diagnostic flag
                if isfield(sensorInfo,'LiGasDiag')
                    irgaGasDiagCol = sensorInfo.LiGasDiag(sensorInfo.LiGasDiag(:,3)==sonHeight,2);
                    gasDiagFlag = output.(tableNames{tble})(:,irgaGasDiagCol);
                    gasDiagFlag(gasDiagFlag > info.diagnosticTest.meanLiGasDiagnosticLimit) = 0;
                    gasDiagFlag(isnan(gasDiagFlag)) = 0;
                else
                    gasDiagFlag = false(size(H2ONanFlag));                    
                end               
                H2OFlag = logical(H2ONanFlag+H2OSpikeFlag+gasDiagFlag);
            else
                irgaH2O = [];
                H2OFlag = [];
            end
            
            if isfield(sensorInfo,'irgaCO2') && ~isempty(sensorInfo.irgaCO2(sensorInfo.irgaCO2(:,3)==sonHeight,2)) % EC150
                irgaCO2col = sensorInfo.irgaCO2(sensorInfo.irgaCO2(:,3)==sonHeight,2);
                irgaCO2SigCol = sensorInfo.irgaCO2sigStrength(sensorInfo.irgaCO2sigStrength(:,3)==sonHeight,2);
                irgaCO2 = data{1,tble}(:,irgaCO2col);  
                CO2sigFlag = output.(tableNames{tble})(:,irgaCO2SigCol);
                CO2sigFlag(CO2sigFlag > info.diagnosticTest.CO2minSignal) = 0;
                CO2sigFlag(isnan(CO2sigFlag)) = 0;
                nanFlagTableName = [tableNames{tble},'NanFlag'];
                spikeFlagTableName = [tableNames{tble},'SpikeFlag'];
                CO2NanFlag = output.(nanFlagTableName)(:,irgaCO2col);
                CO2SpikeFlag = output.(spikeFlagTableName)(:,irgaCO2col);
                CO2Flag = logical(CO2NanFlag+CO2SpikeFlag+gasDiagFlag+CO2sigFlag);
            elseif isfield(sensorInfo,'LiCO2') && ~isempty(sensorInfo.LiCO2(sensorInfo.LiCO2(:,3)==sonHeight,2))%Li7500
                irgaCO2col = sensorInfo.LiCO2(sensorInfo.LiCO2(:,3)==sonHeight,2);
                irgaCO2 = data{1,tble}(:,irgaCO2col)*44; % Mult by 44 to convert from mmol/m^3 to mg/m^3
                nanFlagTableName = [tableNames{tble},'NanFlag'];
                spikeFlagTableName = [tableNames{tble},'SpikeFlag'];
                CO2NanFlag = output.(nanFlagTableName)(:,irgaCO2col);
                CO2SpikeFlag = output.(spikeFlagTableName)(:,irgaCO2col);
                CO2Flag = logical(CO2NanFlag+CO2SpikeFlag+gasDiagFlag);
            else
                irgaCO2 = [];
                CO2Flag = [];
            end
            if isfield(sensorInfo,'KH2O')
                KH2Ocol = sensorInfo.KH2O(sensorInfo.KH2O(:,3)==sonHeight,2);
                KH2O = data{1,tble}(:,KH2Ocol);
                nanFlagTableName = [tableNames{tble},'NanFlag'];
                spikeFlagTableName = [tableNames{tble},'SpikeFlag'];
                H2ONanFlag = output.(nanFlagTableName)(:,KH2Ocol);
                H2OSpikeFlag = output.(spikeFlagTableName)(:,KH2Ocol);
                H2OFlag = logical(H2ONanFlag+H2OSpikeFlag);
            elseif isfield(sensorInfo,'irgaH2O') || isfield(sensorInfo,'LiH2O') % don't delete H2O flag if irga exists!
                KH2O = [];
            else
                KH2O = [];
                H2OFlag = [];
            end
            
            %--------------------- ITERATE THROUGH ALL TIME STEPS
            for jj = 1:N
                if ii == 1
                    % place time stamp in column 1
                    H(jj,1) = t(bp(jj+1));
                    Hlat(jj,1) = t(bp(jj+1));
                    tau(jj,1) = t(bp(jj+1));
                    tke(jj,1) = t(bp(jj+1));
                    LHflux(jj,1) = t(bp(jj+1));
                    CO2flux(jj,1) = t(bp(jj+1));
                    derivedT(jj,1) = t(bp(jj+1));
                    L(jj,1) = t(bp(jj+1));
                    sigma(jj,1) = t(bp(jj+1));
                    if jj == 1
                        Hheader{1} = 'time';
                        HlatHeader{1} = 'time';
                        tauHeader{1} = 'time';
                        tkeHeader{1} = 'time';
                        derivedTheader{1} = 'time';
                        LHfluxHeader{1} = 'time';
                        CO2fluxHeader{1} = 'time';
                        Lheader{1} = 'time';
                        sigmaHeader{1} = 'time';
                    end
                end
                
                % place rho in column 2 of H
                Hheader{2} = 'rho';
                H(jj,2) = rhoAvg(jj);
                
                % place Cp in column 3 of H
                Hheader{3} = 'cp';
                H(jj,3) = 1004.67*(1+0.84*qRefavg(jj)); % pg. 640 of Stull
                
                % find unrotated sonic pertubations
                uP = nandetrend(u(bp(jj)+1:bp(jj+1)),info.detrendingFormat);
                vP = nandetrend(v(bp(jj)+1:bp(jj+1)),info.detrendingFormat);
                wP = nandetrend(w(bp(jj)+1:bp(jj+1)),info.detrendingFormat);
                TsonP = nandetrend(Tson(bp(jj)+1:bp(jj+1)),info.detrendingFormat);
                
                % find rotated sonic pertubations
                uPF_P = nandetrend(uPF(bp(jj)+1:bp(jj+1)),info.detrendingFormat);
                vPF_P = nandetrend(vPF(bp(jj)+1:bp(jj+1)),info.detrendingFormat);
                wPF_P = nandetrend(wPF(bp(jj)+1:bp(jj+1)),info.detrendingFormat);
                
                % store rotated data
                if info.saveRawConditionedData
                    raw.u(bp(jj)+1:bp(jj+1),ii) = u(bp(jj)+1:bp(jj+1));
                    raw.v(bp(jj)+1:bp(jj+1),ii) = v(bp(jj)+1:bp(jj+1));
                    raw.w(bp(jj)+1:bp(jj+1),ii) = w(bp(jj)+1:bp(jj+1));
                    raw.uPF(bp(jj)+1:bp(jj+1),ii) = uPF(bp(jj)+1:bp(jj+1));
                    raw.vPF(bp(jj)+1:bp(jj+1),ii) = vPF(bp(jj)+1:bp(jj+1));
                    raw.wPF(bp(jj)+1:bp(jj+1),ii) = wPF(bp(jj)+1:bp(jj+1));
                    raw.uPF_Prime(bp(jj)+1:bp(jj+1),ii) = uPF_P;
                    raw.vPF_Prime(bp(jj)+1:bp(jj+1),ii) = vPF_P;
                    raw.wPF_Prime(bp(jj)+1:bp(jj+1),ii) = wPF_P;
                    raw.z(ii) = sonHeight;
                end
                
                % find standard deviations of wind vector, sonic and finewire temperature
                if ~isempty(fw) % check for fw first
                    numSigmaVariables = 9;
                    
                    sigma(jj,10+(ii-1)*numSigmaVariables) = nanstd(fw(bp(jj)+1:bp(jj+1))); sigmaHeader{10+(ii-1)*numSigmaVariables} = strcat(num2str(sonHeight),'m :sigma_TFW');
                    if fwFlag(jj); sigma(jj,10+(ii-1)*numSigmaVariables) = nan; end;
                    
                else
                    numSigmaVariables = 8;
                end
                
                sigma(jj,2+(ii-1)*numSigmaVariables) = nanstd(u(bp(jj)+1:bp(jj+1))); sigmaHeader{2+(ii-1)*numSigmaVariables} = strcat(num2str(sonHeight),'m :sigma_u');
                if rotatedSonFlag(jj); sigma(jj,2+(ii-1)*numSigmaVariables) = nan; end;
                
                sigma(jj,3+(ii-1)*numSigmaVariables) = nanstd(v(bp(jj)+1:bp(jj+1))); sigmaHeader{3+(ii-1)*numSigmaVariables} = strcat(num2str(sonHeight),'m :sigma_v');
                if rotatedSonFlag(jj); sigma(jj,3+(ii-1)*numSigmaVariables) = nan; end;
                
                sigma(jj,4+(ii-1)*numSigmaVariables) = nanstd(w(bp(jj)+1:bp(jj+1))); sigmaHeader{4+(ii-1)*numSigmaVariables} = strcat(num2str(sonHeight),'m :sigma_w');
                if unrotatedSonFlag(jj); sigma(jj,4+(ii-1)*numSigmaVariables) = nan; end;
                
                sigma(jj,5+(ii-1)*numSigmaVariables) = nanstd(uPF(bp(jj)+1:bp(jj+1))); sigmaHeader{5+(ii-1)*numSigmaVariables} = strcat(num2str(sonHeight),'m :sigma_uPF');
                if rotatedSonFlag(jj); sigma(jj,5+(ii-1)*numSigmaVariables) = nan; end;
                
                sigma(jj,6+(ii-1)*numSigmaVariables) = nanstd(vPF(bp(jj)+1:bp(jj+1))); sigmaHeader{6+(ii-1)*numSigmaVariables} = strcat(num2str(sonHeight),'m :sigma_vPF');
                if rotatedSonFlag(jj); sigma(jj,6+(ii-1)*numSigmaVariables) = nan; end;
                
                sigma(jj,7+(ii-1)*numSigmaVariables) = nanstd(wPF(bp(jj)+1:bp(jj+1))); sigmaHeader{7+(ii-1)*numSigmaVariables} = strcat(num2str(sonHeight),'m :sigma_wPF');
                if rotatedSonFlag(jj); sigma(jj,7+(ii-1)*numSigmaVariables) = nan; end;
                
                sigma(jj,8+(ii-1)*numSigmaVariables) = nanstd(Tson(bp(jj)+1:bp(jj+1))); sigmaHeader{8+(ii-1)*numSigmaVariables} = strcat(num2str(sonHeight),'m :sigma_Tson');
                if TsonFlag(jj); sigma(jj,8+(ii-1)*numSigmaVariables) = nan; end;
                
                sigma(jj,9+(ii-1)*numSigmaVariables) = nanmean(wPF_P.*TsonP.^2); sigmaHeader{9+(ii-1)*numSigmaVariables} = strcat(num2str(sonHeight),'m :wPFP_TsonP_TsonP');
                if TsonFlag(jj) || rotatedSonFlag(jj); sigma(jj,9+(ii-1)*numSigmaVariables) = nan; end;
                
                % find rotated and unrated momentum flux
                tau(jj,2+(ii-1)*3) = sqrt(nanmean(uP.*wP)^2+nanmean(vP.*wP)^2); tauHeader{2+(ii-1)*3} = strcat(num2str(sonHeight),'m :sqrt(u''w''^2+v''w''^2)');
                if rotatedSonFlag(jj); tau(jj,2+(ii-1)*3) = nan; end;
                
                tau(jj,3+(ii-1)*3) = sqrt(nanmean(uPF_P.*wPF_P)^2+nanmean(vPF_P.*wPF_P)^2); tauHeader{3+(ii-1)*3} = strcat(num2str(sonHeight),'m :sqrt(uPF''wPF''^2+vPF''wPF''^2)');
                if rotatedSonFlag(jj); tau(jj,3+(ii-1)*3) = nan; end;
                
                tau(jj,4+(ii-1)*3) = nanmean(uPF_P.*wPF_P); tauHeader{4+(ii-1)*3} = strcat(num2str(sonHeight),'m :uPF''wPF''');
                if rotatedSonFlag(jj); tau(jj,4+(ii-1)*3) = nan; end;
                
                % find tke
                tke(jj,1+ii) = 1/2*(nanmean(uP.^2)+nanmean(vP.^2)+nanmean(wP.^2)); tkeHeader{1+ii} = strcat(num2str(sonHeight),'m :0.5(u''^2+v''^2+w''^2)');
                if rotatedSonFlag(jj); tke(jj,1+ii) = nan; end;
                
                % find Ts'w' and Ts'wPF' from sonic
                
                H(jj,4+(ii-1)*12) = nanmean(wP.*TsonP); Hheader{4+(ii-1)*12} = strcat(num2str(sonHeight),'m son:Ts''w''');
                if unrotatedSonFlag(jj)||TsonFlag(jj); H(jj,4+(ii-1)*12) = nan; end;
                
                H(jj,5+(ii-1)*12) = nanmean(wPF_P.*TsonP);Hheader{5+(ii-1)*12} = strcat(num2str(sonHeight),'m son:Ts''wPF''');
                if rotatedSonFlag(jj)||TsonFlag(jj); H(jj,5+(ii-1)*12) = nan; end;
                
                Hlat(jj,2+(ii-1)*12) = nanmean(uP.*TsonP);HlatHeader{2+(ii-1)*12} = strcat(num2str(sonHeight),'m son:Ts''u''');
                if rotatedSonFlag(jj)||TsonFlag(jj); Hlat(jj,2+(ii-1)*12) = nan; end;
                
                Hlat(jj,3+(ii-1)*12) = nanmean(uPF_P.*TsonP);HlatHeader{3+(ii-1)*12} = strcat(num2str(sonHeight),'m son:Ts''uPF''');
                if rotatedSonFlag(jj)||TsonFlag(jj); Hlat(jj,3+(ii-1)*12) = nan; end;
                
                Hlat(jj,4+(ii-1)*12) = nanmean(vP.*TsonP);HlatHeader{4+(ii-1)*12} = strcat(num2str(sonHeight),'m son:Ts''v''');
                if rotatedSonFlag(jj)||TsonFlag(jj); Hlat(jj,4+(ii-1)*12) = nan; end;
                
                Hlat(jj,5+(ii-1)*12) = nanmean(vPF_P.*TsonP);HlatHeader{5+(ii-1)*12} = strcat(num2str(sonHeight),'m son:Ts''vPF''');
                if rotatedSonFlag(jj)||TsonFlag(jj); Hlat(jj,5+(ii-1)*12) = nan; end;
                
                % find obukhov length
                kappa = 0.4;
                g = 9.81;
                % T0_L = nanmedian(Tson(bp(j)+1:bp(j+1))) + 273.15; % Tson in K
                T0_L = Tref_Kavg(jj);
                uStarCubed = tau(jj,3+(ii-1)*3).^(3/2); % sqrt(uPF'*wPF')
                H0_L = H(jj,5+(ii-1)*12); % wPF'.*Tson'
                
                L(jj,2+(ii-1)) = -uStarCubed/(kappa*g/T0_L*H0_L); Lheader{2+(ii-1)} = strcat(num2str(sonHeight),'m L:sqrt(uPF''wPF'')^3/2*T_S/(k*g*wPF''Ts'')');
                if rotatedSonFlag(jj)||TsonFlag(jj); L(jj,2+(ii-1)) = nan; end;
                
                % store thetaSonP in raw structure
                if info.saveRawConditionedData
                    raw.sonTs(bp(jj)+1:bp(jj+1),ii) = Tson(bp(jj)+1:bp(jj+1));
                    raw.sonTsPrime(bp(jj)+1:bp(jj+1),ii) = TsonP;
                end
                
                % find Th'w' and Th'wPF' from sonic
                if ~isempty(thetaSonAir)
                    thetaSonAirP = nandetrend(thetaSonAir(bp(jj)+1:bp(jj+1)),info.detrendingFormat);
                    
                    H(jj,6+(ii-1)*12) = nanmean(wP.*thetaSonAirP); Hheader{6+(ii-1)*12} = strcat(num2str(sonHeight),'m son:Th''w''');
                    if unrotatedSonFlag(jj)||TsonFlag(jj); H(jj,6+(ii-1)*12) = nan; end;
                    
                    H(jj,7+(ii-1)*12) = nanmean(wPF_P.*thetaSonAirP);Hheader{7+(ii-1)*12} = strcat(num2str(sonHeight),'m son:Th''wPF''');
                    if rotatedSonFlag(jj)||TsonFlag(jj); H(jj,7+(ii-1)*12) = nan; end;
                end
                
                % find T'w' and T'wPF' from fw
                if ~isempty(fw)
                    fwP = nandetrend(fw(bp(jj)+1:bp(jj+1)),info.detrendingFormat);
                    
                    H(jj,8+(ii-1)*12) = nanmean(wP.*fwP); Hheader{8+(ii-1)*12} = strcat(num2str(sonHeight),'m fw:T''w''');
                    if unrotatedSonFlag(jj)||fwFlag(jj); H(jj,8+(ii-1)*12) = nan; end;
                    
                    H(jj,9+(ii-1)*12) = nanmean(wPF_P.*fwP);Hheader{9+(ii-1)*12} = strcat(num2str(sonHeight),'m fw:T''wPF''');
                    if rotatedSonFlag(jj)||fwFlag(jj); H(jj,9+(ii-1)*12) = nan; end;
                end
                
                % find Th_s'w' and Th_s'wPF' from sonics
                if ~isempty(thetaSon)
                    thetaSonP = nandetrend(thetaSon(bp(jj)+1:bp(jj+1)),info.detrendingFormat);
                    
                    H(jj,10+(ii-1)*12) = nanmean(wP.*thetaSonP); Hheader{10+(ii-1)*12} = strcat(num2str(sonHeight),'m son:Th_s''w''');
                    if unrotatedSonFlag(jj)||TsonFlag(jj); H(jj,10+(ii-1)*12) = nan; end;
                    
                    H(jj,11+(ii-1)*12) = nanmean(wPF_P.*thetaSonP);Hheader{11+(ii-1)*12} = strcat(num2str(sonHeight),'m son:Th_s''wPF''');
                    if rotatedSonFlag(jj)||TsonFlag(jj); H(jj,11+(ii-1)*12) = nan; end;
                    
                    Hlat(jj,6+(ii-1)*12) = nanmean(uP.*thetaSonP); HlatHeader{6+(ii-1)*12} = strcat(num2str(sonHeight),'m son:Th_s''u''');
                    if rotatedSonFlag(jj)||TsonFlag(jj); Hlat(jj,6+(ii-1)*12) = nan; end;
                    
                    Hlat(jj,7+(ii-1)*12) = nanmean(uPF_P.*thetaSonP);HlatHeader{7+(ii-1)*12} = strcat(num2str(sonHeight),'m son:Th_s''uPF''');
                    if rotatedSonFlag(jj)||TsonFlag(jj); Hlat(jj,7+(ii-1)*12) = nan; end;
                    
                    Hlat(jj,8+(ii-1)*12) = nanmean(vP.*thetaSonP); HlatHeader{8+(ii-1)*12} = strcat(num2str(sonHeight),'m son:Th_s''v''');
                    if rotatedSonFlag(jj)||TsonFlag(jj); Hlat(jj,8+(ii-1)*12) = nan; end;
                    
                    Hlat(jj,9+(ii-1)*12) = nanmean(vPF_P.*thetaSonP);HlatHeader{9+(ii-1)*12} = strcat(num2str(sonHeight),'m son:Th_s''vPF''');
                    if rotatedSonFlag(jj)||TsonFlag(jj); Hlat(jj,9+(ii-1)*12) = nan; end;
                end
                
                % find Th'w' and Th'wPF' from fw
                if ~isempty(thetaFw)
                    thetaFwP = nandetrend(thetaFw(bp(jj)+1:bp(jj+1)),info.detrendingFormat);
                    
                    H(jj,12+(ii-1)*12) = nanmean(wP.*thetaFwP); Hheader{12+(ii-1)*12} = strcat(num2str(sonHeight),'m fw:Th''w''');
                    if unrotatedSonFlag(jj)||fwFlag(jj); H(jj,12+(ii-1)*12) = nan; end;
                    
                    H(jj,13+(ii-1)*12) = nanmean(wPF_P.*thetaFwP);Hheader{13+(ii-1)*12} = strcat(num2str(sonHeight),'m fw:Th''wPF''');
                    if rotatedSonFlag(jj)||fwFlag(jj); H(jj,13+(ii-1)*12) = nan; end;
                    
                    Hlat(jj,10+(ii-1)*12) = nanmean(uP.*thetaFwP); HlatHeader{10+(ii-1)*12} = strcat(num2str(sonHeight),'m fw:Th''u''');
                    if rotatedSonFlag(jj)||fwFlag(jj); Hlat(jj,10+(ii-1)*12) = nan; end;
                    
                    Hlat(jj,11+(ii-1)*12) = nanmean(uPF_P.*thetaFwP);HlatHeader{11+(ii-1)*12} = strcat(num2str(sonHeight),'m fw:Th''uPF''');
                    if rotatedSonFlag(jj)||fwFlag(jj); Hlat(jj,11+(ii-1)*12) = nan; end;
                    
                    Hlat(jj,12+(ii-1)*12) = nanmean(vP.*thetaFwP); HlatHeader{12+(ii-1)*12} = strcat(num2str(sonHeight),'m fw:Th''v''');
                    if rotatedSonFlag(jj)||fwFlag(jj); Hlat(jj,12+(ii-1)*12) = nan; end;
                    
                    Hlat(jj,13+(ii-1)*12) = nanmean(vPF_P.*thetaFwP);HlatHeader{13+(ii-1)*12} = strcat(num2str(sonHeight),'m fw:Th''vPF''');
                    if rotatedSonFlag(jj)||fwFlag(jj); Hlat(jj,13+(ii-1)*12) = nan; end;
                    
                    % store fwthetaP in raw structure
                    if info.saveRawConditionedData
                        raw.fwThPrime(bp(jj)+1:bp(jj+1),ii) = thetaFwP;
                        raw.fwTh(bp(jj)+1:bp(jj+1),ii) = thetaFw(bp(jj)+1:bp(jj+1));
                        raw.fwT(bp(jj)+1:bp(jj+1),ii) = fw(bp(jj)+1:bp(jj+1));
                    end
                end
                
                % find VTh'w' and VTh'wPF' from fw
                if ~isempty(VthetaFw)
                    VthetaFwP = nandetrend(VthetaFw(bp(jj)+1:bp(jj+1)),info.detrendingFormat);
                    
                    H(jj,14+(ii-1)*12) = nanmean(wP.*VthetaFwP); Hheader{14+(ii-1)*12} = strcat(num2str(sonHeight),'m fw:VTh''w''');
                    if unrotatedSonFlag(jj)||fwFlag(jj); H(jj,14+(ii-1)*12) = nan; end;
                    
                    H(jj,15+(ii-1)*12) = nanmean(wPF_P.*VthetaFwP);Hheader{15+(ii-1)*12} = strcat(num2str(sonHeight),'m fw:VTh''wPF''');
                    if rotatedSonFlag(jj)||fwFlag(jj); H(jj,15+(ii-1)*12) = nan; end;
                end
                
                % find latent heat flux if it exists at height
                if ~isempty(irgaH2O) || ~isempty(KH2O)
                    
                    % find H2O measurement
                    if ~isempty(irgaH2O)
                        rho_v = irgaH2O; % g/m^3
                        H2Otype = 0;
                        rho_CO2 = irgaCO2; % mg/m^3
                        if isfield(sensorInfo,'irgaH2O') 
                            H2OsensorNumber = find(abs(sensorInfo.irgaH2O(:,3)-sonHeight)<0.2);  % find IRGA number to store raw pertubations
                        end
                        if isfield(sensorInfo,'LiH2O') %&& isempty(H2OsensorNumber)
                            H2OsensorNumber = find(abs(sensorInfo.LiH2O(:,3)-sonHeight)<0.2);  % find IRGA number to store raw pertubations
                        end
                        CO2sensorNumber = H2OsensorNumber;
                    else
                        rho_v = KH2O; % g/m^3
                        H2Otype = 1;
                        H2OsensorNumber = find(abs(sensorInfo.KH2O(:,3)-sonHeight)<0.2);  % find KH2O number to store raw pertubations
                    end
                    
                    % declare constants
                    Mv = 18.0153;   % Molar Mass of H2O (g/mol)
                    Md = 28.97;     % Molar Mass of Dry Air (g/mol)
                    
                    % find heat flux for WPL (wPF_P*TsonP)
                    kinSenFlux = H(jj,5+(ii-1)*12); % K m/s
                    
                    % find the latent heat of vaporizatoin
                    Lv = (2.501-0.00237*(Tref_Kavg(jj)-273.15))*10^3; % Latent Heat of Vaporization (J/g) Stoll P. 641
                    
                    % find H2O pertubations in g/m^3
                    H2Op = nandetrend(rho_v(bp(jj)+1:bp(jj+1)),info.detrendingFormat); %(g/m^3)
                    
                    % store H2OP in raw structure
                    if info.saveRawConditionedData
                        raw.rhovPrime(bp(jj)+1:bp(jj+1),H2OsensorNumber) = H2Op; % (g/m^3)
                        raw.rhov(bp(jj)+1:bp(jj+1),H2OsensorNumber) = rho_v(bp(jj)+1:bp(jj+1)); % (g/m^3)
                    end
                    
                    E = nanmean(wP.*H2Op);     % [g /m^2/s] find evaporation flux from unrotated data
                    EPF = nanmean(wPF_P.*H2Op);  % [g /m^2/s] find evaporation flux from rotated data
                    
                    LHflux(jj,2+(ii-1)*7) = Lv; LHfluxHeader{2+(ii-1)*7} = strcat(num2str(sonHeight),'m Lv(J/g)');
                    
                    LHflux(jj,3+(ii-1)*7) = E; LHfluxHeader{3+(ii-1)*7} = strcat(num2str(sonHeight),'m w'':E(g/m^2s)');
                    if unrotatedSonFlag(jj)||H2OFlag(jj); LHflux(jj,3+(ii-1)*7) = nan; end;
                    
                    LHflux(jj,4+(ii-1)*7) = EPF; LHfluxHeader{4+(ii-1)*7} = strcat(num2str(sonHeight),'m wPF'':E(g/m^2s)');
                    if rotatedSonFlag(jj)||H2OFlag(jj); LHflux(jj,4+(ii-1)*7) = nan; end;
                    
                    % apply O2 correction to KH2O
                    if H2Otype % see: https://www.eol.ucar.edu/content/corrections-sensible-and-latent-heat-flux-measurements 
                               % EVERYTHING ELSE IS DEFINITELY WRONG!!! THIS LOOKS RIGHT
                        ko = -0.0045; % effective O2 absorption coefficient.  Value from: Tanner et al., 1993
                        kw = -0.153;  % effective H2O absorption coefficient.  Value from instrument calibration.  -0.153 is from 2012 calibration of the Pardyjak lab KH2O
                        CkO = 0.23*ko/kw;
                        O2correction = CkO*rhodAvg(jj)/Tref_Kavg(jj)*kinSenFlux*1000;  % multiply by 1000 to put in g/m^2/s
                        E = E+O2correction;
                        EPF = EPF+O2correction;
                        LHflux(jj,7+(ii-1)*7) = Lv*E; LHfluxHeader{7+(ii-1)*7} = strcat(num2str(sonHeight),'m O2 no WPL,w'' (W/m^2)');
                        LHflux(jj,8+(ii-1)*7) = Lv*EPF; LHfluxHeader{8+(ii-1)*7} = strcat(num2str(sonHeight),'m O2 no WPL,wPF'' (W/m^2)');
                        
                        % create WPL headers to include O2 correction
                        LHfluxHeader{5+(ii-1)*7} = strcat(num2str(sonHeight),'m WPL and O2, w'' (W/m^2)');
                        LHfluxHeader{6+(ii-1)*7} = strcat(num2str(sonHeight),'m WPL and O2, wPF'' (W/m^2)');
                    end
                    
                    % apply WPL Corrections (E.C. by Marc Aubinet 97).  Headers are created above for KH2O and below for EC150
                    LHflux(jj,5+(ii-1)*7) = 1000*Lv*(1+Md/Mv*rhovAvg(jj)/rhodAvg(jj))*(E./1000+(rhovAvg(jj)/Tref_Kavg(jj))*kinSenFlux); 
                    if unrotatedSonFlag(jj)||H2OFlag(jj); LHflux(jj,5+(ii-1)*7) = nan; end;
                    
                    LHflux(jj,6+(ii-1)*7) = 1000*Lv*(1+Md/Mv*rhovAvg(jj)/rhodAvg(jj))*(EPF./1000+(rhovAvg(jj)/Tref_Kavg(jj))*kinSenFlux); 
                    if unrotatedSonFlag(jj)||H2OFlag(jj); LHflux(jj,6+(ii-1)*7) = nan; end;
                    
                    
                    if ~H2Otype  % find CO2 flux and create H2O WPL headers
                        
                        % create WPL, H2O headers 
                        LHfluxHeader{5+(ii-1)*7} = strcat(num2str(sonHeight),'m WPL, w'' (W/m^2)');
                        LHfluxHeader{6+(ii-1)*7} = strcat(num2str(sonHeight),'m WPL, wPF'' (W/m^2)');
                        
                        rho_CO2p = nandetrend(rho_CO2(bp(jj)+1:bp(jj+1)),info.detrendingFormat)./1e6; % kg/m^3
                        rho_CO2avg = nanmean(rho_CO2(bp(jj)+1:bp(jj+1)))/1e6; % kg/m^3
                        evapFlux = LHflux(jj,6+(ii-1)*7)/Lv/1000; % WPL, wPF'' (kg/m^2s)
                        
                        CO2flux(jj,2+(ii-1)*4) = nanmean(wP.*rho_CO2p); CO2fluxHeader{2+(ii-1)*4} = strcat(num2str(sonHeight),'m w'':CO2(kg/m^2s)');
                        if unrotatedSonFlag(jj)||CO2Flag(jj); CO2flux(jj,2+(ii-1)*4) = nan; end;
                        
                        CO2flux(jj,3+(ii-1)*4) = nanmean(wPF_P.*rho_CO2p); CO2fluxHeader{3+(ii-1)*4} = strcat(num2str(sonHeight),'m wPF'':CO2(kg/m^2s)');
                        if rotatedSonFlag(jj)||CO2Flag(jj); CO2flux(jj,3+(ii-1)*4) = nan; end;
                        
                        CO2flux(jj,4+(ii-1)*4) = CO2flux(jj,2+(ii-1)*4) + Md/Mv*(rho_CO2avg/rhodAvg(jj))*evapFlux + (1+Md/Mv*rhovAvg(jj)/rhodAvg(jj))*(rho_CO2avg/Tref_Kavg(jj))*kinSenFlux; CO2fluxHeader{4+(ii-1)*4} = strcat(num2str(sonHeight),'m WPL,w'':CO2(kg/m^2s)');
                        
                        CO2flux(jj,5+(ii-1)*4) = CO2flux(jj,3+(ii-1)*4) + Md/Mv*(rho_CO2avg/rhodAvg(jj))*evapFlux + (1+Md/Mv*rhovAvg(jj)/rhodAvg(jj))*(rho_CO2avg/Tref_Kavg(jj))*kinSenFlux; CO2fluxHeader{5+(ii-1)*4} = strcat(num2str(sonHeight),'m WPL,wPF'':CO2(mol/m^2s)');
                        % store CO2P in raw structure
                        if info.saveRawConditionedData
                            raw.rhoCO2Prime(bp(jj)+1:bp(jj+1),CO2sensorNumber) = rho_CO2p;
                            raw.rhoCO2(bp(jj)+1:bp(jj+1),CO2sensorNumber) = rho_CO2(bp(jj)+1:bp(jj+1));
                        end
                    end
                end
                
                
            end
        catch err
            message = strcat(err.message,'@ line',num2str(err.stack.line),' Problem with sonic at ',num2str(sonHeight),' m will be skipped');
            warning(message)
            if isempty(output.warnings{1})
                output.warnings{1,1} = message;
            else
                output.warnings{end+1,1} = message;
            end
        end
    end
    %------------- STORE OUTPUTS
    flag = logical(any(H,1)+isnan(H(1,:)));
    output.H = H(:,flag);
    output.Hheader = Hheader(flag);
    
    flag = logical(any(Hlat,1)+isnan(Hlat(1,:)));
    output.Hlat = Hlat(:,flag);
    output.HlatHeader = HlatHeader(flag);
    
    flag = logical(any(tau,1)+isnan(tau(1,:)));
    output.tau = tau(:,flag);
    output.tauHeader = tauHeader(flag);
    
    flag = logical(any(tke,1)+isnan(tke(1,:)));
    output.tke = tke(:,flag);
    output.tkeHeader = tkeHeader(flag);
    
    flag = logical(any(sigma,1)+isnan(sigma(1,:)));
    output.sigma = sigma(:,flag);
    output.sigmaHeader = sigmaHeader(flag);
    
    flag = logical(any(L,1)+isnan(L(1,:)));
    output.L = L(:,flag);
    output.Lheader = Lheader(flag);
    
    if size(derivedT,2) > 1
        flag = logical(any(derivedT,1)+isnan(derivedT(1,:)));
        output.derivedT = derivedT(:,flag);
        output.derivedTheader = derivedTheader(flag);
    end
    
    if size(LHflux,2) > 1
        flag = logical(any(LHflux,1)+isnan(LHflux(1,:)));
        output.LHflux = LHflux(:,flag);
        output.LHfluxHeader = LHfluxHeader(flag);
    end
    
    if size(CO2flux,2) > 1
        flag = logical(any(CO2flux,1)+isnan(CO2flux(1,:)));
        output.CO2flux = CO2flux(:,flag);
        output.CO2fluxHeader = CO2fluxHeader(flag);
    end
catch err
    message = strcat(err.message,'@ line ',num2str(err.stack.line),' UNABLE TO FIND FLUXES AT All HEIGHTS');
    warning(message)
    raw = [];
    if isempty(output.warnings{1})
        output.warnings{1,1} = message;
    else
        output.warnings{end+1,1} = message;
    end
end
end