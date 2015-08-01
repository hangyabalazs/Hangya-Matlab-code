function somsfcoh(cellids)
%SOMSFCOH   Spike-field coherence.
%   SOMSFCOH calculates spike-field coherence between local field potential
%   (LFP) and single cell spiking data. The spike train (point process) is
%   convolved with a Gaussian window (49 ms, SD = 6 ms) using FILTFILT
%   (no shifting of the peaks occur). Signals are downsampled at 200 Hz.
%   Magnitude squared coherence is calculated using MSCOHERE with a
%   64-points Hanning window. The resulting coherence spectrum is plotted
%   and saved.
%
%   SOMSFCOH(I) uses cells defined by the index set I.
%   SOMSFCOH(CELLIDS) uses cells defined by the list of CellIDs (CELLIDS).
%
%   See also MSCOHERE and SOMSTA.

% Pass the control to the user in case of error
dbstop if error

% Directories
global DATAPATH
fs = filesep;
resdir = [DATAPATH 'SOM' fs 'sfcoh' fs];  % results directory
fnm = 'SFCOH_matrices.mat';   % filename for saving the result matrices

% Input argument check
if nargin < 1
    loadcb   % load CellBase
    cellids = CELLIDLIST; 
else
    if isnumeric(cellids)
        loadcb   % load CellBase
        cellids = CELLIDLIST(cellids);
    end
end
if ischar(cellids)   % convert 'cellids' to cell
    cellids = {cellids};
end

% STA loop
NumCells = length(cellids);   % number of cells
for iC = 1:NumCells
    cellid = cellids{iC};
    [Cxy fc Pxx fpx Pyy fpy Hc Hp] = sfcohmain(cellid);   % calculate spike-field coherence in a loop
    
    % Save
    cellidt = regexprep(cellid,'\.','_');
    ff1 = [resdir 'SFCOH1_' cellidt '.fig'];   % save coherence figure
    saveas(Hc,ff1)
    ff2 = [resdir 'SFCOH2_' cellidt '.fig'];   % save coherence figure with individual spectra
    saveas(Hp,ff2)
    ff3 = [resdir 'SFCOH_' cellidt '.mat'];   % save coherence matrices
    save(ff3,'Cxy','fc','Pxx','fpx','Pyy','fpy')
    close all
end

% -------------------------------------------------------------------------
function [Cxy fc Pxx fpx Pyy fpy Hc Hp] = sfcohmain(cellid)

% LFP full path
[r s tetrodename] = cellid2tags(cellid);   %#ok<*ASGLU> % get tetrode number
matname = cellid2fnames(cellid,'cont',tetrodename);   % .mat filename for LFP
[pathname filename extension] = fileparts(matname);   % parse filename
cscname = fullfile(pathname,[filename '.ncs']);   % filename for the Neuralynx CSC file

% Load LFP
[TimeStamp, ChanNum, SampleFrequency, NumValSamples, Samples, NlxHeader] = ...
    Nlx2MatCSC(cscname,[1 1 1 1 1],1,1,1);  % mPFC LFP
lfp_orig = -1 * Samples(:);   % invert Neuralynx data
sr_orig = 512 / mean(diff(TimeStamp)) * 10^6;    % originaml sampling rate (changed later)
dt_orig = 1 / sr_orig;   % original timestep
starttime = TimeStamp(1) / 10^6;   % start time of recording
time_orig0 = (0:length(lfp_orig)-1) * dt_orig + starttime;   % generate time vector
time_orig = repmat(TimeStamp/10^6,512,1) + repmat((0:511)'*dt_orig,1,length(TimeStamp));
time_orig = time_orig(:)';

% Load Events
fn = [pathname '\EVENTS.mat'];
load(fn);         % load converted Neuralynx events

% Use the last recording from file
start_recording = find(Events_EventIDs==0);   % find recording start time(s)
if length(start_recording) > 1   % more than one starts
    disp([cellid ': Multiple recording sessions found. Longest recording is used.'])
    endinx = [start_recording(2:end)-1 length(Events_EventIDs)];   % last event of each recording
    lenrec = Events_TimeStamps(endinx) - Events_TimeStamps(start_recording);  % length of the recording
    chinx = lenrec==max(lenrec);   % index for the longest recording
    pst = Events_TimeStamps(start_recording(chinx));   % use longest recording
    starttime = time_orig(find(time_orig>=pst,1,'first'));   % start time of last recording
    linx = time_orig > pst;  % indices of last recording
    lfp_orig = lfp_orig(linx);   % restrict LFP
    time_orig = time_orig(linx);   % restrict time vector
end

% Resample
sr = 1000;  % new sampling rate
dt = 1 / sr;   % new time step
[p q] = rat(sr/sr_orig);    % resample LFP at 1000 Hz
lfp = resample(lfp_orig,p,q);
time = (0:length(lfp)-1) * dt + starttime;   % generate new time vector

% Load spike train
limit_spikes = [0 100000];   % include max 100000 spikes
try
    tseg = findSegs3(cellid,'segfilter','stim_excl');  % find time segments
    unit = extractSegSpikes(cellid,tseg);   % find spikes in the time segments
catch ME
    disp([cellid ': findSegs3 exited with the following error.'])
    disp(ME.message)
    disp([cellid ': All spikes are used.'])
    unit = loadcb(cellid,'SPIKES');   % use all spikes
end
if any(unit>time(end)) | any(unit<time(1))   %#ok<OR2> % drop spikes out of LFP time range
    noso = sum(unit>time(end)|unit<time(1));   % number of out-of-range spikes
    unit(unit>time(end)|unit<time(1)) = [];   % restrict to LFP time range
    disp([cellid ': ' num2str(noso) ' spikes out of LFP time range.'])
end
if length(unit) > limit_spikes(2);      % crop if too long to avoid out of memory
    unit = unit(1:limit_spikes(2));
end
NumSpikes = length(unit);
vdisc = round((unit-starttime)*sr) + 1;

% Spike Density Function
sdf = gconv(vdisc,length(lfp));  % convolution with a Gaussian kernel

% Downsample
dsr = 100;   % downsample at 100 Hz (alternative: downsample at 200 Hz, 256-point Hanning-window)
cnst = sr / dsr;   % use integer ratios and downsample without interpolation
sdf2 = sdf(1:cnst:end);   % downsampled spike density function
lfp2 = lfp(1:cnst:end);   % downsampled LFP

% Linear detrending
dsdf = detrend(sdf2);
dlfp = detrend(lfp2);

% Coherence
[Cxy fc] = mscohere(dsdf,dlfp,hanning(64),[],[],dsr);   % magnitude squared coherence

% Periodogram
[Pxx fpx] = periodogram(dsdf,[],[],dsr);
[Pyy fpy] = periodogram(dlfp,[],[],dsr);

% Plot
plinx = fc < 50;  % plot below 50 Hz
Hc = figure;
plot(fc(plinx),Cxy(plinx))   % freq. vs. coherence
xlabel('Frequency [Hz]')
ylabel('Coherence')
title('Spike-field coherence')

Hp = figure;
subplot(223)   % unit spectrum
plinx = fpx < 50;  % plot below 50 Hz
plot(fpx(plinx),Pxx(plinx))   % freq. vs. spectral power
xlabel('Frequency [Hz]')
ylabel('Power spectral density')
title('Unit spectrum')

subplot(222)   % LFP spectrum
plinx = fpy < 50;  % plot below 50 Hz
plot(fpy(plinx),Pyy(plinx))   % freq. vs. spectral power
xlabel('Frequency [Hz]')
ylabel('Power spectral density')
title('LFP spectrum')

subplot(224)   % coherence
plinx = fc < 50;  % plot below 50 Hz
plot(fc(plinx),Cxy(plinx))   % freq. vs. coherence
xlabel('Frequency [Hz]')
ylabel('Coherence')
title('Spike-field coherence')

% -------------------------------------------------------------------------
function sdf = gconv(vdisc,lenlfp)

% Window function
gwn = gausswin(49,4);   % 49 ms Gaussian window; SD = ((n-1) / 2) / alpha = 6 ms;
gwn = gwn / sum(gwn);

% Convolution
psunit = zeros(1,lenlfp);
psunit(vdisc) = 1;
sdf = filtfilt(gwn,1,psunit);