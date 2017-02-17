% Utah Turbulence in Environmental Studies Process and Analysis Code (UTESpac) Created by: Derek Jensen,
% derek591@gmail.com Version 4.0 Version Date: 1 June 2016
% 
% About:
% 
% UTESpac is designed specifically for use with Campbell Scientific dataloggers and accompanying LoggerNet software with
% native support for
% 
% 	Sonic Anemometers:  RMYOUNG 8100, Campbell Sci CSAT3 Open Path Gas Analyzers:  Licor 7500, Campbell Sci EC150 and
% 	IRGASON, Krypton Hygrometers Finewire thermocouples for heat flux computations Propeller Anemometers Mean
% 	meteorological sensors (e.g. T/RH, Pressure, Solar, cup anemometers, etc.)
% 	
% UTESpac expects 24 or 48-hr CSV tables, quality controls the data and then computes means, fluxes, variances and
% derived temperatures (potential temperature, virtual potential temperature) and stores the output in a MatLab
% structure or NetCDF file.
% 
% Steps for Use:
% 
% 1.  Convert Campbell Binary files to csv files using the Card Convert Program in LoggerNet
% 	Options: File Processing - Use Time, set to 2 days 00 h under Time Settings
% 			 File naming - Use TimeDate Filenames and Append to Last File if multiple site files exist Array CSV Options
% 			 Timestamp Options - Include year, day, hour/minutes, seconds, don't include midnight is 2400, Array ID,
% 			 Array Datalogger Format = Hour/Minutes and Seconds
%  
%  2.  Create a folder for the individual site.  The folder name needs to be preceded by the keyword "site".  E.g. for a
%  site named Playa the folder name is sitePlaya.  Place the .csv files within the site folder
%  
%  3.  Create a subfolder named output, this is where the output data will be stored
%  
%  4.  Create header files for each data table The syntax is ‘tableName’_header.dat (e.g. "Playa_1HZ_header.dat",
%  "Playa_20HZ_header.dat).  Note that ‘tableName’ must be consistent with the .csv tableNames created in step 1.  The
%  header file is a single line .dat, comma delimited file containing variable names and heights for all columns within
%  the respective data table.  The header file is 3 columns shorter than the .csv data file.  This is because UTESpac
%  immediately calculates the serial date numbers from the date vectors (columns 1 – 4) contained in the data tables.
%  The serial dates are stored in column 1 and columns 2 – 4 are deleted, thus becoming consistent with the header file.
%  The easiest way to create the header file is with Card Convert.  Create an ASCII T0A5 file, there is no need to run
%  the whole binary file, simply stop the conversion immediately and only a few hundred lines will be created.  Open the
%  file in a text editor and delete all lines outside of the variable headers (typically line 5).  The variable names
%  within the header and the sensor templates (defined on lines 155-169) must be consistent.  The template is used by
%  UTESpac to identify specific sensors in the header.  The rules for creating the template and header variable names
%  are:
% 	- The template and variable name are the exact same except the sensor height is replaced with the wildcard '*' in
% 	the template.  e.g. template = 'Ux_*', header variable name = 'Ux_0.5', 'Ux_10'
%    - The sensor height must be the last numeric value in the header variable name - All sensors (with exception of
%    solar and battery) need an associated height in meters - Heights within the header variable name at a given tower
%    height need to exactly match. e.g. 'FW_5','Ux_5','RH_5'
%  
%  5.  If a global planar fit is used, a PFinfo structure, containing global planar fit coefficients, will be stored in
%  the site folder.  There is no need to do anything with it.  Note: For the Global Planar Fit, there must be 1 and only
%  1 set of 5 minute, local planar fit data.  That is, the global planar fit will fail if there is
%  '5minAvg_LPF_linDetrend' and '5minAvg_LPF_constDetrend' in the output folder.  There must be one or the other (it
%  doesn't matter which!).
% 
% 6.  Fill out the information section of the code (lines 56 - 116) and run the code.  A full example study is included
% in UTESpac.zip
% 
% 7.  Use getData(), structFill() and structConcat() to produce complete (no missing days) datasets over the full
% experiment.  See example
%
%Eric's comments questions:
%1. Need to have a siteInfo.m where the site specific information (e.g. orientation of sonic
% is added) 
%2.Need showcell with UTESPac


%% ----------------------------------------------USER INFORMATION-------------------------------------------------------

% current UTESpac Version
info.UTESpacVersion = '4.0';

% enter root folder where site* folders are located
info.rootFolder = 'C:\Users\derek\desktop\';

% enter averaging period in minutes.  Must yield an integer when dividied into 60 (e.g. 1, 2, 5, 10, 20, 30)
info.avgPer = 5;

% save QC'd raw data tables in .mat format
info.saveRawConditionedData = true;

% save netCDF file
info.saveNetCDF = false;

% save .csv files
info.saveCSV = false;

% enter detrending format ('constant' or 'linear')
info.detrendingFormat = 'linear';

% select 'local' or 'global' planar fit, 'local' computes coeffiecients from local file only, 'global' computes
% user-defined, multi-sector, multi-datebin coefficients from all site data - the sector and datebins are defined
% graphically when the code is executed - for 'global' calculations, all data must first be run with a 'local' planar
% fit and 5-min averaging
info.PF.globalCalculation = 'local';

% recalulate global PF coefficients if 'global' calculation is used
info.PF.recalculateGlobalCoefficients = false;

% select averaging period for global PF calculation, if used - local PF calculation runs with average specified in
% info.avgPer
info.PF.avgPer = 30;

% chose max and min wind magintudes for global PF calculation, if used - local PF calculation uses all velocity data
info.PF.globalCalcMaxWind = 12;
info.PF.globalCalcMinWind = 0.5;

% give reference specific humidity if humidity measurement does not exist
info.qRef = 12; % (g/kg)
% ---

% --- DATA CONDITIONING OPTIONS See: Vickers and Mahrt 97 and envsupport.licor.com/help/EddyPro3/Content/Topics/Despiking_Raw_Stat_Screening.htm
% 1. spike removal and spike test
info.spikeTest.maxRuns = 20;  % max number of time to run through despike algorithm
info.spikeTest.windowSizeFraction = 1;  % fraction that is multiplied by info.avgPer to get width of averaging window.  Must be <= 1
info.spikeTest.maxConsecutiveOutliers = 10;  % maximum number of consecutive spikes that will be removed.  Longer segments are considered physical
info.spikeTest.maxPercent = 2; % max acceptable percentage of spikes per averaging period
% number of std's to define a spike
info.spikeTest.spikeDef.u = 3.5;
info.spikeTest.spikeDef.v = 3.5;
info.spikeTest.spikeDef.w = 5;
info.spikeTest.spikeDef.Tson = 3.5;
info.spikeTest.spikeDef.fw = 3.5;
info.spikeTest.spikeDef.irgaCO2 = 3.5;
info.spikeTest.spikeDef.irgaH2O = 3.5;
info.spikeTest.spikeDef.KH2O = 3.5;
info.spikeTest.spikeDef.cup = 3.5;
info.spikeTest.spikeDef.birdSpd = 3.5;
info.spikeTest.spikeDef.otherInstrument = 5;

% 2. absolute limits - NaN values outside of limits
info.absoluteLimitsTest.u = [-50, 50]; % m/s
info.absoluteLimitsTest.v = [-50, 50];
info.absoluteLimitsTest.w = [-10, 10];
info.absoluteLimitsTest.Tson = [-20, 80];  % deg C
info.absoluteLimitsTest.fw = [-20, 80];
info.absoluteLimitsTest.irgaCO2 = [0, 1500]; % mg/m^3
info.absoluteLimitsTest.irgaH2O = [0, 50];  % g/m^3
info.absoluteLimitsTest.KH2O = [0, 50];
info.absoluteLimitsTest.cup = [0 50];
info.absoluteLimitsTest.birdSpd = [0 50];

% 3. wind direction
% enter +/- envelope around around tower orientation to define 'bad wind direction'
info.windDirectionTest.envelopeSize = 20;

% 4. excessive NaN test
info.nanTest.maxPercent = 55;

% 5. diagnostic tests
info.diagnosticTest.H2OminSignal = 0.7;
info.diagnosticTest.CO2minSignal = 0.7;
info.diagnosticTest.meanGasDiagnosticLimit = 0.1;
info.diagnosticTest.meanSonicDiagnosticLimit = 50;
info.diagnosticTest.meanLiGasDiagnosticLimit = 220;  % Full strength is 255, less indicates problems.  See manual. 
% ---

% specify sensor templates that must match those found in the header files, place wild card (*) over sensor height [expected units]
template.u = 'Ux_*'; % sonic u  --   [m/s]
template.v = 'Uy_*'; % sonic v  --   [m/s]
template.w = 'Uz_*'; % sonic w  --   [m/s]
template.Tson = 'T_Sonic_*'; % sonic T  --   [C or K]
template.sonDiagnostic = 'diagnostic_*'; % sonic diagnostic  --  [-]
template.fw = 'FW_*'; % sonic finewires to be used for Eddy Covariance  --  [C]
template.RH = 'RH_*'; % slow response relative humidity for virtual temperature calculation  --  [Fract or %]
template.T = 'Temp_*'; % slow response temperature  --  [C]
template.P = 'Pressure_*'; % pressure  --  [kPa or mBar]
template.irgaH2O = 'H2O_*'; % for use with Campbell EC150 and IRGASON.  WPL corrections applied  --  [g/m^3]
template.irgaH2OsigStrength = 'H2Osig_*'; % EC150 Signal Strength  --  [-]
template.irgaCO2 = 'CO2_*'; % for use with Campbell EC150 and IRGASON.  WPL corrections applied  --  [mg/m^3]
template.irgaCO2sigStrength = 'CO2sig_*'; % EC150 Signal Strength  --  [-]
template.irgaGasDiag = 'gas_diag_*'; % EC150 gas (CO2 and H2O) diagnostic, 0-> Okay  --  [-]
template.LiH2O = 'LiH2O_*'; % for use with Licor 7500.  WPL corrections applied  --  [mmol/mol]
template.LiCO2 = 'LiCO2_*'; % for use with Licor 7500.  WPL corrections applied  --  [mmol/mol]
template.LiGasDiag = 'Li_gas_diag_*'; % Li7500 gas (CO2 and H2O) diagnostic >~230 -> Okay  --  [-]
template.KH2O = 'KH2O_H2O_*'; % for use with KH2Os.  WPL and O2 corrections applied  --  [g/m^3]
template.cup = 'cup_*';  % wind speed from cup anemometers  --  [m/s]
template.birdSpd = 'wbSpd_*';  % wind speed from prop anemometer  --  [m/s]
template.birdDir = 'wbDir_*';  % wind direction from vain or prop anemometer  --  [deg]
%% ------------------------------------------------ANALYSIS-------------------------------------------------------------
% find data files and load headers
[headers, dataFiles, tableNames, info] = findFiles(info);

% find instrument locations
[sensorInfo, info] = findInstruments(headers,template,info);

% find global PF coefficients
if strcmp(info.PF.globalCalculation,'global')
    PFinfo = findGlobalPF(info,template,sensorInfo);
else
    PFinfo = [];
end

for i = 1:size(dataFiles,1)
    try
        
        % load files
        [data, dataInfo] = loadData(dataFiles(i,:),i,size(dataFiles,1),info,tableNames);
        
        % store serial date in column 1 of data tables and delete columns 2-4
        [data, dataInfo, info] = findSerialDate(data, dataInfo, info);
        
        % condition data
        [data, output] = conditionData(data,info,tableNames,template,headers);
        
        % store raw averages in output structure
        output =  avg(data,info,tableNames,output,headers,sensorInfo);
        
        % store wind statistics
        output = windStats(output,sensorInfo,tableNames,info);

        % apply local or global (info.PFrange option) planar fit and yaw rotations(Wilczak et. al, 2000)
        [rotatedSonicData, output, dataInfo] = sonicRotation(output,data,sensorInfo,info,dataInfo,tableNames,PFinfo);

        % find sensible heat flux, momentum flux, latent heat flux, CO2 flux
        [output, raw] = fluxes(data, rotatedSonicData, info, output, sensorInfo,tableNames);

        % save data in the root folder
        output = saveData(info,output,dataInfo,headers, tableNames, raw, template);
                
    catch err
        
        warning('Problem with %s @ line %g:%s',dataFiles{i,1},err.stack(2,1).line, err.message)
        
    end
end
gong();