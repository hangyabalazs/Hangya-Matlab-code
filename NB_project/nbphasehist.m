function nbphasehist(cellids)
%NBPHASEHIST   Phase histogram.
%   NBPHASEHIST calculates theta phase histogram. Stimulation segments are
%   excluded from the analysis. The resulting phase histogram plot is saved
%   along with all phase values and descriptive phase statistics.
%
%   Phase is calculated via Hilbert-transform of the LFP filtered between 4
%   and 12 Hz. Filtering is performed by FIR filter applying bidirectional
%   zero-phase filtering with a filter order of 4096. Cycles with a length
%   out of the frequency range or with an amplitude lower than mean + 2 *
%   SD are discarded.
%
%   NBPHASEHIST(I) uses cells defined by the index set I.
%   NBPHASEHIST(CELLIDS) uses cells defined by the list of CellIDs
%   (CELLIDS).
%
%   See also PHASEHIST.

% Pass the control to the user in case of error
dbstop if error

% Directories
global DATAPATH
fs = filesep;
resdir = [DATAPATH 'NB' fs 'phasehist_' fs];  % results directory
fnm = 'PHASE_matrices.mat';   % filename for saving the result matrices

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

% Phase histogram loop
NumCells = length(cellids);   % number of cells
for iC = 1:NumCells
    cellid = cellids{iC};
    [hang hmvl ftm angs p_rayleigh p_rao H] = phasemain(cellid);   % calculate phase distribution in a loop
    
    % Save
    cellidt = regexprep(cellid,'\.','_');
    ff1 = [resdir 'PHASE_' cellidt '.fig'];   % save phase histogram figure
    saveas(H,ff1)
    ff2 = [resdir 'PHASE_' cellidt '.mat'];   % save phase distribution matrices
    save(ff2,'hang','hmvl','ftm','p_rayleigh','p_rao','angs')
    close all
end

% -------------------------------------------------------------------------
function [hang hmvl ftm angs p_rayleigh p_rao H] = phasemain(cellid)

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

% Use the longest recording from file
start_recording = find(strcmp(Events_EventStrings,'Starting Recording'));   % find recording start time(s)
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
    tseg = findSegs3(cellid,'segfilter','stim_excl_nb');  % find time segments
    unit = extractSegSpikes(cellid,tseg);   % find spikes in the time segments
catch ME
    disp('findSegs3 exited with the following error.')
    disp(ME.message)
    disp('All spikes are used.')
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

% Phase histogram
[hang hmvl ftm angs p_rayleigh p_rao H] = somphase(lfp,vdisc,sr,4,12);   % call phase histogram code [0.5-4 for delta]

% -------------------------------------------------------------------------
function [hang, hmvl, ftm, bahee, p_rayleigh, p_rao, H] = somphase(eeg,vdisc,sr,fr1,fr2)
%SOMPHASE   Phase histogram.
%   [A R FTM ANGS P1 P2 H] = SOMPHASE(EEG,VD,SR,FR1,FR2) calculates phase
%   histograms for neuronal spikes (VD) relative to EEG sampled at SR.
%   Phase is calculated via Hilbert-transform of the EEG filtered between
%   FR1 and FR2 Hz. Output arguments are mean phase (A), mean vector length
%   (R), first trigonometric moment (FTM), all phase values (ANGS), p
%   values of significant phase locking with Rayleight-test (P1) and Rao's
%   spacing test (P2) and handle of the histogram figure (H).
%
%   See also CZPHASE2.

% Input argument check
vdisc = vdisc(:)';   % convert to row vector

% Filtering EEG
nqf = sr / 2;
flt = fir1(4096,[fr1 fr2]/nqf,'band');      % bandpass filtering
feeg = filtfilt(flt,1,eeg);
feeg = (feeg - mean(feeg)) / std(feeg);

% Hilbert transformation of the EEG
ahee = angle(hilbert(feeg));
aheedeg = ahee * (180 / pi);

% Check criteria:
% 1. discard cicles with EEG amp. lower then 2SD
% 2. discard cicles shorter then 1/fr2 s
% 3. discard cicles longer then 1/fr1 s
% fn = find(-diff(ahee)>2*pi-0.3);
fn0 = valuecrossing(1:length(ahee),ahee',0,'up');
fn = round(fn0);
sd = std(feeg);
inx = [];
vdisc(vdisc<0|vdisc>length(eeg)) = [];
for k = 1:length(fn)-1
    seeg = feeg(fn(k)+1:fn(k+1));
    axs = max(seeg) - min(seeg);
    if (axs < 2 * sd)  || (fn(k+1) - fn(k) > 1 / fr1 * sr) || (fn(k+1) - fn(k) < 1 / fr2 * sr)
        inx = [inx find(vdisc>fn0(k)&vdisc<fn0(k+1))];
    end
end
vdisc(inx) = [];

% Phase angles - Hilbert
bahee = ahee(round(vdisc));
n = length(bahee);
ftm = sum(exp(1).^(i*bahee)) / n;    % first trigonometric moment
hang = angle(ftm);   % mean angle
hmvl = abs(ftm);     % mean resultant length

% Plot phase histogram
angs = bahee;
edges = -180:20:180;     % edges for phase histogram
cnts = (edges(1:end-1) + edges(2:end)) / 2;
[nm,xout] = histc(angs*180/pi,edges);   % phase histogram
nm = nm(1:end-1);
H = figure;
B = bar(cnts,nm'/length(angs));
set(B,'FaceColor',[0.16 0.38 0.27])
y_lim = ylim;
axis([-200 200 y_lim(1) y_lim(2)])
str = ['\it{Mean angle: }' '\bf ' num2str(hang*180/pi)];
text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','black')
str = ['\it{Mean resultant length: }' '\bf ' num2str(hmvl)];
text(60,y_lim(2)-1.5*(y_lim(2)-y_lim(1))/6,str,'Color','black')
str = ['\it{n: }' '\bf ' num2str(length(angs))];
text(60,y_lim(2)-2*(y_lim(2)-y_lim(1))/6,str,'Color','black')
if ~isempty(bahee)      % significance testing
    [Z,p_rayleigh,U,p_rao] = b_rao(bahee(:)');
    if p_rayleigh < 0.001
        clr_ray = 'red';
    else
        clr_ray = 'black';
    end
    if p_rao(2) <= 0.05
        clr_rao = 'red';
    else
        clr_rao = 'black';
    end
    str = ['\it{Raylegh''s p = }' '\bf ' num2str(p_rayleigh)];
    text(60,y_lim(2)-2.5*(y_lim(2)-y_lim(1))/6,str,'Color',clr_ray)
    str = ['\it{Rao''s p < }' '\bf ' num2str(p_rao(2))];
    text(60,y_lim(2)-3*(y_lim(2)-y_lim(1))/6,str,'Color',clr_rao)
end