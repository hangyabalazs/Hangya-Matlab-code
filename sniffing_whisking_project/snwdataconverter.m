function snwdataconverter(rrat,rdate,csc1,csc2,msgn)
%SNWDATACONVERTER   Convert Neuralynx data to Matlab.
%   SNWDATACONVERTER(RAT,DATE,CSC1,CSC2,MSGN) accepts input arguments
%   determining the animal (RAT) and session (DATE) information, which data
%   channels are used (CSC1 and CSC2 for sniffing and whisking,
%   respectively) and a scalar determining whether sniffing is inverted
%   (MSGN=-1) or not (MSGN=1). It saves the neural data, the position data
%   and Neuralynx events in .mat format.
%
%   See also SNWSNIFFPLOT2.

% Input arguments
error(nargchk(0,5,nargin))
if nargin < 1
    rrat = 'R7';    % animal
end
if nargin < 2
    rdate = '2008-09-03_15-59-08';   % session
end
if nargin < 3
    csc1 = 'CSC4.ncs';  % respiration (sniffing) channel
end
if nargin < 4
    csc2 = 'CSC6.ncs';  % EMG (whisking) channel
end
if nargin < 5
    msgn = -1;   % control the polarity of the sniffing signal
end
dsr = 1000;   % downsample rate

% Define results directory
global DATADIR
resdir = [DATADIR 'HSW\mat2\'];

% Read data
fn = ([DATADIR 'HSW\' rrat '\' rdate '\']);
[Sniffing.TimeStamp, Sniffing.ChanNum, Sniffing.SampleFrequency, ...
    Sniffing.NumValSamples, Sniffing.Samples, Sniffing.NlxHeader] = ...
    Nlx2MatCSC([fn csc1],[1 1 1 1 1],1,1,1);  % sniffing
sniff = msgn * Sniffing.Samples(:);
sr_sniff = 512 / (Sniffing.TimeStamp(2) - Sniffing.TimeStamp(1)) * 10^6;    % originaml sampling rate (changed later)
dt_sniff = 1 / sr_sniff;
starttime1 = Sniffing.TimeStamp(1) / 10^6;
time_sniff = (0:length(sniff)-1) * dt_sniff + starttime1;
[Whisking.TimeStamp, Whisking.ChanNum, Whisking.SampleFrequency, ...
    Whisking.NumValSamples, Whisking.Samples, Whisking.NlxHeader] = ...
    Nlx2MatCSC([fn csc2],[1 1 1 1 1],1,1,1);  % whisking
whisk = Whisking.Samples(:);
sr_whisk = 512 / (Whisking.TimeStamp(2) - Whisking.TimeStamp(1)) * 10^6;    % originaml sampling rate for EMG (different from that for respiration in R7)
dt_whisk = 1 / sr_whisk;
starttime2 = Whisking.TimeStamp(1) / 10^6;
time_whisk = (0:length(whisk)-1) * dt_whisk + starttime2;
timediff = starttime2 - starttime1;
% sniff = linterp(time_sniff,sniff,time_whisk);

% Resample
[p q] = rat(dsr/sr_sniff);    % resample respiration at 1000 Hz
sniff = resample(sniff,p,q);
time_sniff = (0:length(sniff)-1) / dsr + starttime1;
[p q] = rat(dsr/sr_whisk);    % resample EMG at 1000 Hz
whisk = resample(whisk,p,q);
time_whisk = (0:length(sniff)-1) / dsr + starttime2;
sr = dsr;  % new sampling rate
if timediff > 0.0001       % deal with misaligned data
    if ~isequal(sr,1000)
        error('Data realignment only works for samplin rate = 1000 Hz.')
    end
    leadtime = timediff * sr;    % lead of sniff in ms
    sniff(1:round(leadtime)) = [];  % drop lead: now they co-start
    time_sniff(1:round(leadtime)) = [];
elseif timediff < -0.0001
    if ~isequal(sr,1000)
        error('Data realignment only works for samplin rate = 1000 Hz.')
    end
    lagtime = (-1) * timediff * sr;    % lag of sniff in ms
    whisk(1:round(lagtime)) = [];  % drop lead of whisk: now they co-start
    time_whisk(1:round(lagtime)) = [];
end
ln = min(length(sniff),length(whisk));
sniff = sniff(1:ln);
time_sniff = time_sniff(1:ln);
whisk = whisk(1:ln);
time_whisk = time_whisk(1:ln);   % now they co-terminate
if abs(max(time_sniff-time_whisk)) > 0.001   % time vectors should agree with 1 ms tolerance
    error('Time vector problem.')
else
    time = time_sniff;
end

% Events
[Events.EventTimeStamps, Events.EventIDs, Events.Nttls, ...
    Events.Extras, Events.EventStrings Events.NlxHeader] = ...
    Nlx2MatEV([fn 'Events.nev'],[1 1 1 1 1],1,1,1);         % load Neuralynx events

% Position data
[Position.TimeStamps, Position.ExtractedX, Position.ExtractedY, ...
    Position.ExtractedAngle, Position.Targets, Position.Points, Position.NlxHeader] =...
    Nlx2MatVT([fn 'VT1.nvt'],[1 1 1 1 1 1],1,1);   % load Neuralynx tracking file
xpos = Position.ExtractedX;
ypos = Position.ExtractedY;
time_pos = Position.TimeStamps;

if any(xpos.*ypos) 
    
    % Position vector correction
    [xpos ypos time_pos] = poscorrection(xpos,ypos,time_pos);
    
    % Resample position vector
    sr_pos = 1e6 / mean(diff(time_pos));    % original sampling rate (changed later)
    dt_pos = 1 / sr_pos;
    starttime_pos = time_pos(1) / 10^6;
    [p q] = rat(dsr/sr_pos);    % resample position vector at 1000 Hz
    xpos = resample(xpos,p,q);
    ypos = resample(ypos,p,q);
    time_pos = starttime_pos:1/dsr:starttime_pos+(length(xpos)-1)/dsr;
    timediff = starttime_pos - time(1);
    if timediff > 0.05       % deal with misaligned data
        if ~isequal(sr,1000)
            error('Data realignment only works for sampling rate = 1000 Hz.')
        end
        leadtime = timediff * sr;    % lead of sniff in ms
        sniff(1:round(leadtime)) = [];  % drop lead: now they co-start
        whisk(1:round(leadtime)) = [];
        time_sniff(1:round(leadtime)) = [];
        time_whisk(1:round(leadtime)) = [];
    elseif timediff < -0.05
        if ~isequal(sr,1000)
            error('Data realignment only works for sampling rate = 1000 Hz.')
        end
        lagtime = (-1) * timediff * sr;    % lag of sniff in ms
        xpos(1:round(lagtime)) = [];  % drop lead of whisk: now they co-start
        ypos(1:round(lagtime)) = [];
        time_pos(1:round(lagtime)) = [];
    end
    ln = min(length(sniff),length(xpos));
    sniff = sniff(1:ln);
    time_sniff = time_sniff(1:ln);
    whisk = whisk(1:ln);
    time_whisk = time_whisk(1:ln);
    xpos = xpos(1:ln);
    ypos = ypos(1:ln);
    time_pos = time_pos(1:ln);   % now they co-terminate
    if abs(max(time_sniff-time_pos)) > 0.2   % time vectors should agree with 100 ms tolerance
        error('Time vector problem.')
    else
        time = time_pos;
    end
end

% Save
fnm = [resdir 'DATA_' rrat '_' rdate '.mat'];   % save variables
save(fnm,'Sniffing','Whisking','Events','Position',...
    'sniff','whisk','time','xpos','ypos','time_pos',...
    'csc1','csc2','msgn','rrat','rdate')

% -------------------------------------------------------------------------
function [X Y TP] = poscorrection(x_pos,y_pos,x_pos_ts)

% Position vector correction
bns = 8;
TP2 = x_pos_ts;
dx = diff(x_pos);       % leave out the outliers
fdx = find(abs(dx)>bns&abs([dx(2:end) 0])>bns&abs(dx-[dx(2:end) 0])>2*bns);
x_pos(fdx+1) = [];
y_pos(fdx+1) = [];
x_pos_ts(fdx+1) = [];
dy = diff(y_pos);
fdy = find(abs(dy)>bns&abs([dy(2:end) 0])>bns&abs(dy-[dy(2:end) 0])>2*bns);
x_pos(fdy+1) = [];
y_pos(fdy+1) = [];
x_pos_ts(fdy+1) = [];

inxs = x_pos > 0 & y_pos > 0;   % leave out (0;0) points
X = x_pos(inxs);
Y = y_pos(inxs);
TP = x_pos_ts(inxs);

% Interpolate the deleted outlier timepoints.
X = interp1(TP,X,TP2);
Y = interp1(TP,Y,TP2);
TP = TP2;