function AVG = simpleAvg(input,avgPer,varargin)
%simpleAvg.m finds the time average of a matrix based on the serial date number. The first averaging period will be from
%start_hour:00 to start_hour:avg_per and will end at start_hour:00 on the last day of data.  It is more generalized than
%avg.m.  !!The matrix should first be filled with timeStampFill since averaging bins are determined with linspace.m
%rather than find.m!! Optional arguments: returnTimestamps - true(default) time stamps returned; false time stamps
%ommited from averaged matrix WDpntr & WSpntr: logical array indicating wind direction and wind speed for vector average
%wind speed calculation

if nargin == 3
    returnTimeStamps = varargin{1};
elseif nargin == 4
    WDpntr = varargin{1};
    WSpntr = varargin{2};
elseif nargin == 5
    returnTimeStamps = varargin{1};
    WDpntr = varargin{2};
    WSpntr = varargin{3};
end

% find serial date column
tCol = find(input(1,:) > datenum(2000,1,1) & input(1,:) < datenum(2020,1,1),1,'last');

% find total number of averaging bins
N = round(ceil((input(end,tCol))-floor(input(1,tCol)))/(avgPer/(24*60)));

% ensure that time stamp spacing is consistent within 0.5 seconds and that the input has more rows than N
if isempty(tCol) || N > size(input,1)  || range(diff(input(:,tCol))) > datenum(0,0,0,0,0,0.5)
    warning('Time Stamp Spacing Inconsistent or Input Matrix Incomplete.  Matrix Was Not Averaged')
    AVG = input;
    return
end

% find averaging breakpoints. size(bp,1) = size(N,1) + 1.  Round() is being used as a precaution
bp = round(linspace(0,size(input,1),N+1));

% initialize averaged matrix
AVG = nan(N,size(input,2));

% iterate through time steps.  bp(i) marks the end of the previous bin
for i = 1:N
    AVG(i,:) = nanmean(input(bp(i)+1:bp(i+1),:),1);
    AVG(i,tCol) = input(bp(i+1),tCol);
    
    % apply averaging for windbirds if WBpntr and WSpntr are input
    if exist('WDpntr','var') && exist('WSpntr','var')        
        WS = input(bp(i)+1:bp(i+1),WSpntr);
        WD = deg2rad(input(bp(i)+1:bp(i+1),WDpntr));
        v = nanmean(WS.*sin(WD),1);
        u = nanmean(WS.*cos(WD),1);
        AVG(i,WDpntr) = mod(atan2(v,u)*180/pi,360);
    end
        
end

if exist('returnTimeStamps','var')
    if ~returnTimeStamps
        AVG = AVG(:,1:end-1);
    end
end