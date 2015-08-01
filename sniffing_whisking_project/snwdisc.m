function snwdisc(rrat,rdate,csc1,csc2,msgn)
%SNWDISC   Discrimination of whisking events.
%   SNWDISC performs threshold based discrimination of whisking events. the
%   user is prompted to choose an appropriate threshold. Raw whisking data
%   is resampled at 1 kHz and high-pass filtered at 10 Hz.
%
%   SNWDISC(RAT,DATE,CSC1,CSC2,MSGN) accepts input arguments determining
%   the animal (RAT) and session (DATE) information, which data channels
%   are used (CSC1 and CSC2 for sniffing and whisking, respectively) and a
%   scalar determining whether sniffing is inverted (MSGN=-1) or not
%   (MSGN=1).
%
%   See also SNWDISC_CALL and UDISC.

% Input arguments
error(nargchk(0,5,nargin))
if nargin < 1
    rrat = 'R3';
end
if nargin < 2
    rdate = '2008-07-09_14-46-47_R3';
end
if nargin < 3
    csc1 = 'CSC4.ncs';
end
if nargin < 4
    csc2 = 'CSC5.ncs';
end
if nargin < 5
    msgn = 1;
end

% Define results directory
dbstop if error
global DATAPATH
resdir = [DATAPATH 'SniffWhisk\disc\'];

% Data import
dsr = 1000;   % downsample rate
plotting = true;
global DATADIR
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
flt = fir1(4096,10/nqf,'high');
fwhisk = filtfilt(flt,1,whisk);     % highpass filter whisking at 10 Hz

% Discriminate whisking
[vdisc_rawwhisk disc_thr] = b_udisc(fwhisk');

% Save
fnm = [resdir 'DISC_' rrat '_' rdate '.mat'];   % save variables
save(fnm,'vdisc_rawwhisk','disc_thr','csc1','csc2','msgn','rrat','rdate')
% keyboard