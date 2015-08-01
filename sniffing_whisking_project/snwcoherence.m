function snwcoherence(rrat,rdate,csc1,csc2,msgn,thr,wl)
%SNWCOHERENCE   Coherence between whisking and respiration cycle.
%   SNWCOHERENCE calculates mean-square coherence between whiking root-mean
%   square and sniffing. Sniffing is low-pass filtered at 30 Hz and
%   whisking is high-pass filtered at 10 Hz.
%
%   SNWCOHERENCE(RAT,DATE,CSC1,CSC2,MSGN,THR,WL) accepts input arguments
%   determining the animal (RAT) and session (DATE) information, which data
%   channels are used (CSC1 and CSC2 for sniffing and whisking,
%   respectively), a scalar determining whether sniffing is inverted
%   (MSGN=-1) or not (MSGN=1), threshold for RMS discrimination (THR) and
%   window size for RMS in ms (WL).
%
%   See also SNWSNIFFPLOT_CALL and MSCOHERE..

% Input arguments
error(nargchk(0,7,nargin))
if nargin < 1
    rrat = 'P5';    % animal
end
if nargin < 2
    rdate = 'P5 2005-11-7';   % session
end
if nargin < 3
    csc1 = 'CSC4.ncs';  % respiration (sniffing) channel
end
if nargin < 4
    csc2 = 'CSC6.ncs';  % EMG (whisking) channel
end
if nargin < 5
    msgn = 1;   % control the polarity of the sniffing signal
end
if nargin < 6
    thr = 85;  % threshold for RMS discrimination
end
if nargin < 7
    wl = 20;    % windows for RMS
end

% Define results directory
global DATAPATH
global DATADIR
resdir = [DATAPATH 'SniffWhisk\coherence3\'];

% Load discriminated raw whisking
fnme = [DATAPATH 'SniffWhisk\disc\DISC_' rrat '_' rdate '.mat'];
load(fnme)
vdisc_rawwhisk = vdisc_rawwhisk / wl;

% Data import
plotting = false;
dsr = 1000;   % downsample rate

fn = ([DATADIR 'HSW\' rrat '\' rdate '\']);
[TimeStamp, ChanNum, SampleFrequency, NumValSamples, Samples, NlxHeader] = ...
    Nlx2MatCSC([fn csc1],[1 1 1 1 1],1,1,1);  % sniffing
sniff = msgn * Samples(:);
sr_sniff = 512 / mean(diff(TimeStamp)) * 10^6;    % original sampling rate (changed later)
dt_sniff = 1 / sr_sniff;
starttime1 = TimeStamp(1) / 10^6;
time_sniff = (0:length(sniff)-1) * dt_sniff + starttime1;
[TimeStamp, ChanNum, SampleFrequency, NumValSamples, Samples, NlxHeader] = ...
    Nlx2MatCSC([fn csc2],[1 1 1 1 1],1,1,1);  % whisking
whisk = Samples(:);
sr_whisk = 512 / mean(diff(TimeStamp)) * 10^6;    % original sampling rate for EMG (different from that for respiration in R7)
dt_whisk = 1 / sr_whisk;
starttime2 = TimeStamp(1) / 10^6;
time_whisk = (0:length(whisk)-1) * dt_whisk + starttime2;
timediff = starttime2 - starttime1;
% sniff = linterp(time_sniff,sniff,time_whisk);

% Resample
[p q] = rat(dsr/sr_sniff);    % resample respiration at 1000 Hz
sniff = resample(sniff,p,q);
[p q] = rat(dsr/sr_whisk);    % resample EMG at 1000 Hz
whisk = resample(whisk,p,q);
sr = dsr;  % new sampling rate
if timediff > 0.0001       % deal with misaligned data
    if ~isequal(sr,1000)
        error('Data realignment only works for samplin rate = 1000 Hz.')
    end
    leadtime = timediff * sr;    % lead of sniff in ms
    sniff(1:round(leadtime)) = [];  % drop lead: now they co-start
elseif timediff < -0.0001
    if ~isequal(sr,1000)
        error('Data realignment only works for samplin rate = 1000 Hz.')
    end
    lagtime = (-1) * timediff * sr;    % lag of sniff in ms
    whisk(1:round(lagtime)) = [];  % drop lead of whisk: now they co-start
end
ln = min(length(sniff),length(whisk));
sniff = sniff(1:ln);
whisk = whisk(1:ln);   % now they co-terminate

% Deal with individual artifacts
if isequal(rdate,'2008-09-03_16-42-14')
    whisk = whisk(1:2.25*10^5);
    sniff = sniff(1:2.25*10^5);
elseif isequal(rdate,'2008-09-21_21-08-00')
    whisk = whisk(1:1.54*10^6);
    sniff = sniff(1:1.54*10^6);
end

% Filter
nqf = sr / 2;
flt = fir1(1024,30/nqf,'low');
fsniff = filtfilt(flt,1,sniff);     % lowpass filter sniffing at 30 Hz
flt = fir1(4096,10/nqf,'high');
fwhisk = filtfilt(flt,1,whisk);     % highpass filter whisking at 10 Hz

% Root mean square of whisking
lenw = length(fwhisk);
lle = floor(lenw/wl) * wl;
wh2 = reshape(fwhisk(1:lle),wl,lle/wl);
rms_whisk = sqrt(sum(wh2.^2)) / sqrt(wl);
rms_whisk = rms_whisk';
sniff2 = fsniff(1:lle);
sniff2 = sniff2(round(wl/2):wl:end);   % downsampled sniff aligned to whisking RMS

% Coherence
% [Cxy f] = mscohere(sniff,whisk,hamming(1024),[],[],dsr);
% [Cxy f] = mscohere(rms_whisk,sniff2,hanning(1024),512,1024,sr/wl);    % in folder coherence2 
[Cxy f] = mscohere(rms_whisk,sniff2,hanning(128),[],[],sr/wl);

% Save
fnm = [resdir 'COHERENCE_' rrat '_' rdate '.mat'];   % save variables
save(fnm,'Cxy','f','wl','csc1','csc2','msgn','rrat','rdate')

% Plot
plinx = f < 50;
H = figure;
plot(f(plinx),Cxy(plinx))

% Save
fnm = [resdir 'COHERENCE_' rrat '_' rdate '.fig'];   % save figure
saveas(H,fnm)