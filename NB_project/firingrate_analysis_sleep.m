function firingrate_analysis_sleep(cellid,varargin)
%FIRINGRATE_ANALYSIS_SLEEP   Firing rate.
%   FR = FIRINGRATE_ANALYSIS_SLEEP(CELLID) calculates firing rate in sleep 
%   sessions. LFP spectrum, delta to all power ratio, temporal evolution of
%   firing rate and the animals speed, speed - firing rate association are
%   plotted.
%
%   Note that the sleep session MAY INCLUDE LIGHT_STIMULATION EPISODES!
%
%   See also FIRINGRATE_ANALYSIS and NBFIRINGRATE.

%   Balazs Hangya, Cold Spring Harbor Laboratory
%   1 Bungtown Road, Cold Spring Harbor
%   balazs.cshl@gmail.com
%   27-Aug-2012

%   Edit log: BH, 8/27/13

% Default arguments
prs = inputParser;
addRequired(prs,'cellid',@iscellid)
parse(prs,cellid,varargin{:});
g = prs.Results;

% Load stimulus events
try
    SE = loadcb(cellid,'StimEvents');   % load events
    SS = loadcb(cellid,'STIMSPIKES');   % load prealigned spikes
catch ME
    disp('There was no stimulus protocol for this session.')
    error(ME.message)
end

% Load absoulte spike times
SP = loadcb(cellid);   % absoulte spike times

% Firing rate binned to 30s segments
seglen = 30;        % 30 sec. long segments
len = SP(end) - SP(1);   % length of the recording
lim1 = SP(1):seglen:SP(end)-seglen;   % segment boundaries
lim2 = lim1 + seglen;
fr_time = (lim1 + lim2) / 2;   % middle time points for FR vector
NumSeg = floor(len/seglen);   % number of segments
fr = nan(1,NumSeg);
for iS = 1:NumSeg
    segspikes = SP(SP>lim1(iS)&SP<lim2(iS));  % spikes in the segment
    fr(iS) = length(segspikes) / seglen;  % FR in the segment
end

% Position data
[animalID sessionID] = cellid2tags(cellid);
fullpth = fullfile(getpref('cellbase','datapath'),animalID,sessionID);
fn = fullfile(fullpth,'VT1.nvt');
[TimeStamps, ExtractedX, ExtractedY, ExtractedAngle, Targets, Points, NlxHeader] =...
    Nlx2MatVT(fn,[1 1 1 1 1 1],1,1);   % load position data
TimeStamps = TimeStamps / 1e6;   % convert time stamps to seconds

% Spatial firing pattern
% [irhst,rhst,rhst2,thst,thst2] = place(SP,ExtractedX,ExtractedY,TimeStamps,...
%     'BinSize',16,'Display',true);   % rate map

% Position correction
[xpos ypos time_pos] = poscorrection(ExtractedX,ExtractedY,TimeStamps);
if ~isequal(time_pos,TimeStamps)  % internal check
    error('Bug.')
end

% Resample position vector
dt_pos = mean(diff(time_pos));
sr_pos = 1 / dt_pos;    % original sampling rate (changed later)

% Speed
rx = diff(prctile(xpos,[0.25 99.75]));   % x range, allowing some outliers
ry = diff(prctile(ypos,[0.25 99.75]));   % y range, allowing some outliers
xpos2 = xpos / rx * 27.4;  % x position in cm
ypos2 = ypos / ry * 16.5;  % y position in cm
spd = speed(xpos2,ypos2,sr_pos);

% [p q] = rat(1000/sr_pos);    % upsample speed vector at 1000 Hz
% spd2 = resample(spd,p,q,0);   % nearest-neighbour interpolation (step function), very fast
time_pos2 = time_pos(1):0.001:time_pos(end);
spd2 = interp1(time_pos(1:end-1),spd,time_pos2);   % linear interpolation for upsampling to 1000 Hz
zunit = zeros(size(time_pos2));   % digitized unit
zunit(round((SP-time_pos2(1))*1000)) = 1;

% Speed-firing rate association
maxsp = ceil(max(spd));
bns = 2;   % speed bin size
bno = floor(maxsp/bns);  % number of speed bins
[spksp timesp frsp] = deal(nan(1,bno));
bndrs = 0:bns:maxsp;   % bin boundaries
bncnt = (bndrs(1:end-1) + bndrs(2:end)) / 2;  % speed bin centers
for k = 1:bno   % loop though speed bins
    inx = spd2 >= bndrs(k) & spd2 < bndrs(k+1);   % index array for the speed bin
    spksp(k) = sum(zunit(inx));   % number of spikes in the speed bin
    timesp(k) = sum(inx) / 1000;   % cumulative time in s
    frsp(k) = spksp(k) / timesp(k);   % FR in the speed bin
end
figure
plot(bncnt,frsp)
xlabel('speed (cm/s)')
ylabel('firing rate (Hz)')

% Speed in half-min bins
sp_tmb = nan(1,NumSeg);
for iS = 1:NumSeg   % speed binned according to time
    segspeed = spd(time_pos>lim1(iS)&time_pos<lim2(iS));  % speed in the segment
    sp_tmb(iS) = mean(segspeed);  % mean speed in the segment
end
figure
plot(fr_time,sp_tmb)
hold on
plot(fr_time,fr,'r')

% Load LFP
cscname = fullfile(fullpth,'CSC5.ncs');   % filename for the Neuralynx CSC file
[TimeStamps_CSC, ~, ~, ~, Samples, NlxHeader] = ...
    Nlx2MatCSC(cscname,[1 1 1 1 1],1,1,1);  % mPFC LFP
TimeStamps_CSC = TimeStamps_CSC / 1e6;   % convert to seconds
starttime = TimeStamps_CSC(1);   % first time stamp
Samples = Samples(1:32:end);   % downsample
lfp = -1 * Samples(:)';   % flip data (Nlx convention is inverted)
clear Samples
sr = 16 / (mean(diff(TimeStamps_CSC)));   % sampling rate
dt = 1 /sr;  % sampling time
% time = TimeStamps_CSC(1):dt:TimeStamps_CSC(1)+dt*(length(lfp)-1);  % time vector (even)
time_orig = repmat(TimeStamps_CSC,16,1) + repmat((0:15)'*dt,1,length(TimeStamps_CSC));
time_orig = time_orig(:)';   % time vector (accurate, reflecting uneven sampling of Digilynx)

% Load Events
fn = [fullpth '\EVENTS.mat'];
load(fn);         % load converted Neuralynx events
figure;
plot(Events_TimeStamps,Events_Nttls)

% Resample LFP
dsr = 200;  % new sampling rate
ddt = 1 / dsr;   % new time step
[p q] = rat(dsr/sr);    % resample LFP at 500 Hz
lfp2 = resample(lfp,p,q);
time2 = (0:length(lfp2)-1) * ddt + starttime;   % generate new time vector
figure;
plot(time2,lfp2)

% Wavelet
[wave_eeg, f] = eegwavelet(lfp2,dsr);        % EEG wavelet
figure
imagesc(wave_eeg)
ff = round(f*100) / 100;
wavetime = round(time2*100) / 100;
b_rescaleaxis('Y',ff)
b_rescaleaxis('X',wavetime)

% Delta power
wave_eeg(wave_eeg>0.1) = 0.1;
ratio = deltaratio(wave_eeg,f);
figure
plot(time2(1:100:end),smooth(ratio(1:100:end),'linear',101))

% Output


keyboard

% -------------------------------------------------------------------------
function [wave,f] = eegwavelet(dat,sr)

% Prepare for wavelet transformation
variance = std(dat) ^ 2;
dat = (dat - mean(dat)) / sqrt(variance) ;
n = length(dat);
dt = 1 / sr;
pad = 0;  % No padding!
dj = 0.02;    
j1 = ceil((1/dj) * log2(n/2));
j1 = ceil(j1);
j = (0:j1);
s0 = 2 * dt; 
s = s0 .* 2 .^ (j * dj);
omega0 = 6;
c = 4 * pi / (omega0 + sqrt(2+omega0^2));
fperiod = c .* s;
f = 1 ./ fperiod;
lag1 = 0.72;
param = -1;
mif = 1;          %minimal intresting frequency
mis = find(f>mif);
mis = mis(end);     %maximal intristing scale
mother = 'Morlet';

% Wavelet transformation
[wave,period,scale,coi] = b_wavelet_pow3(dat,dt,pad,dj,s0,j1,mother,param,mis);

% --------------------------------------------------------------------------------
function ratio = deltaratio(power,f)

% Computing delta power
fnd = find(f>3);  % upper freq. boundary
pwind_delta1 = fnd(end);
fnd = find(f<1.2);   % lower freq. boundary
pwind_delta2 = fnd(1);
deltapower = power(pwind_delta1:pwind_delta2,:);   % delta band
sumdeltapower = sum(deltapower);   % delta power
clear deltapower

% Compute ratio
ratio = sumdeltapower ./ sum(power);   % ratio of delta power to all power

% -------------------------------------------------------------------------
function [X Y TP] = poscorrection(x_pos,y_pos,x_pos_ts)

% Position vector correction
bns = 8;   % spatial bin size
TP2 = x_pos_ts;   % for later interpolation
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

% -------------------------------------------------------------------------
function spd = speed(xpos,ypos,sr)

% Speed
scxp = smooth(xpos,'linear',51);   % smoothed x positions
scyp = smooth(ypos,'linear',51);   % smoothed y positions
dx = diff(scxp);   % derivative of x positions
dy = diff(scyp);   % derivative of y positions
spd = sqrt(dx.^2+dy.^2) * sr;   % speed