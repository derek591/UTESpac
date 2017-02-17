Utah Turbulence in Environmental Studies Process and Analysis Code (UTESpac) 
Created by: Derek Jensen and Eric Pardyjak
derek591@gmail.com 
Version 4.1 
Version Date: 15 January 2017

About:

UTESpac is designed specifically for use with Campbell Scientific dataloggers and accompanying LoggerNet software with
native support for 

      Sonic Anemometers:  RMYOUNG 8100, Campbell Sci CSAT3 Open Path Gas Analyzers:  Licor 7500, Campbell Sci EC150 and 
      IRGASON, Krypton Hygrometers Finewire thermocouples for heat flux computations Propeller Anemometers Mean
      meteorological sensors (e.g. T/RH, Pressure, Solar, cup anemometers, etc.)
   
UTESpac expects 24 or 48-hr CSV tables, quality controls the data and then computes means, fluxes, variances and 
derived temperatures (potential temperature, virtual potential temperature) and stores the output in a MatLab
structure or NetCDF file.

Steps for Use:

1.  Convert Campbell Binary files to csv files using the Card Convert Program in LoggerNet
      Options: File Processing - Use Time, set to 2 days 00 h under Time Settings
                       File naming - Use TimeDate Filenames and Append to Last File if multiple site files exist Array CSV Options
                       Timestamp Options - Include year, day, hour/minutes, seconds, don't include midnight is 2400, Array ID, 
                       Array Datalogger Format = Hour/Minutes and Seconds
 
 2.  Create a folder for the individual site.  The folder name needs to be preceded by the keyword "site".  E.g. for a
 site named Playa the folder name is sitePlaya.  Place the .csv files within the site folder
 
 3.  Create a subfolder named output, this is where the output data will be stored
 
 4.  Create header files for each data table The syntax is <91>tableName<92>_header.dat (e.g. "Playa_1HZ_header.dat",
 "Playa_20HZ_header.dat).  Note that <91>tableName<92> must be consistent with the .csv tableNames created in step 1.  The 
 header file is a single line .dat, comma delimited file containing variable names and heights for all columns within
 the respective data table.  The header file is 3 columns shorter than the .csv data file.  This is because UTESpac
 immediately calculates the serial date numbers from the date vectors (columns 1 <96> 4) contained in the data tables.
 The serial dates are stored in column 1 and columns 2 <96> 4 are deleted, thus becoming consistent with the header file.
 The easiest way to create the header file is with Card Convert.  Create an ASCII T0A5 file, there is no need to run 
 the whole binary file, simply stop the conversion immediately and only a few hundred lines will be created.  Open the 
 file in a text editor and delete all lines outside of the variable headers (typically line 5).  The variable names
 within the header and the sensor templates (defined on lines 155-169) must be consistent.  The template is used by
 UTESpac to identify specific sensors in the header.  The rules for creating the template and header variable names
 are:
   - The template and variable name are the exact same except the sensor height is replaced with the wildcard '*' in
   the template.  e.g. template = 'Ux_*', header variable name = 'Ux_0.5', 'Ux_10'
   - The sensor height must be the last numeric value in the header variable name - All sensors (with exception of
   solar and battery) need an associated height in meters - Heights within the header variable name at a given tower
   height need to exactly match. e.g. 'FW_5','Ux_5','RH_5'
 
 5.  If a global planar fit is used, a PFinfo structure, containing global planar fit coefficients, will be stored in
 the site folder.  There is no need to do anything with it.  Note: For the Global Planar Fit, there must be 1 and only
 1 set of 5 minute, local planar fit data.  That is, the global planar fit will fail if there is
 '5minAvg_LPF_linDetrend' and '5minAvg_LPF_constDetrend' in the output folder.  There must be one or the other (it 
 doesn't matter which!).

6.  Fill out the information section of the code (lines 56 - 116) and run the code.  A full example study is included
in UTESpac.zip

7.  Use getData(), structFill() and structConcat() to produce complete (no missing days) datasets over the full
experiment.  See example
