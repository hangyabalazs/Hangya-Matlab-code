function nblfp(cellids)

%SOMSTA   Spike Triggered Average.
%   SOMSTA calculates Spike Triggered Average in one second windows.
%   Stimulation segments are excluded from the analysis. STA is visualized
%   along with 95% confidence limits. Confidence interval is calculated by
%   a bootstrap algorithm using frequency-adjusted random Poisson units.
%
%   SOMSTA(I) uses cells defined by the index set I.
%   SOMSTA(CELLIDS) uses cells defined by the list of CellIDs (CELLIDS).
%
%   See also SOMPVCCGCALL2.

% Pass the control to the user in case of error
dbstop if error

% Directories
global DATAPATH
fs = filesep;
resdir = [DATAPATH 'NB' fs 'sta' fs];  % results directory
fnm = 'STA_matrices.mat';   % filename for saving the result matrices

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
    [st sta stats H] = stamain(cellid);   % calculate STA in a loop
    H2 = gcf;
    
    % Save
    cellidt = regexprep(cellid,'\.','_');
    ff1 = [resdir 'STA_' cellidt '.fig'];   % save STA figure
    saveas(H,ff1)
    ff2 = [resdir 'STAimage_' cellidt '.jpg'];   % save STA figure
    saveas(H2,ff2)
    ff3 = [resdir 'STA_' cellidt '.mat'];   % save STA matrices
    save(ff3,'sta','stats')
    close all
end

% -------------------------------------------------------------------------
function [st sta stats H] = stamain(cellid)

% LFP full path
[r s tetrodename] = cellid2tags(cellid);   %#ok<*ASGLU> % get tetrode number
matname = cellid2fnames(cellid,'cont',tetrodename);   % .mat filename for LFP
[pathname filename extension] = fileparts(matname);   % parse filename
cscname = fullfile(pathname,[filename '.ncs']);   % filename for the Neuralynx CSC file

% Load LFP
[TimeStamp, ChanNum, SampleFrequency, NumValSamples, Samples, NlxHeader] = ...
    Nlx2MatCSC(cscname,[1 1 1 1 1],1,1,1);  % mPFC LFP
lfp_orig = -1 * Samples(:);   % invert Neuralynx data
sr_orig = 512 / (mean(diff(TimeStamp))) * 10^6;    % originaml sampling rate (changed later)
dt_orig = 1 / sr_orig;   % original timestep
starttime = TimeStamp(1) / 10^6;   % start time of recording

% Load Events
fn = [pathname '\EVENTS.mat'];
load(fn);         % load converted Neuralynx events

% Use the longest recording from file
start_recording = find(strcmp(Events_EventStrings,'Starting Recording'));   % find recording start time(s)
if length(start_recording) > 1   % more than one starts
    disp([cellid ': Multiple recording sessions found. Longest recording is used.'])
    time_orig = repmat(TimeStamp/10^6,512,1) + repmat((0:511)'*dt_orig,1,length(TimeStamp));
    time_orig = time_orig(:)';
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
clear lfp_orig
time = (0:length(lfp)-1) * dt + starttime;   % generate new time vector

% Determine time window
wn = 1 * sr;    % 1 s window

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
TE = loadcb(cellid,'TrialEvents');
SE = loadcb(cellid,'StimEvents');
stimon = TE.TrialStart(TE.Hit==1|TE.Miss==1) + TE.StimulusOn(TE.Hit==1|TE.Miss==1);
stimon = TE.TrialStart(TE.Hit==1|TE.FalseAlarm==1) + TE.StimulusOn(TE.Hit==1|TE.FalseAlarm==1);
stimon = TE.TrialStart(TE.Hit==1) + TE.StimulusOn(TE.Hit==1);
stimon = TE.TrialStart(TE.StimulusDuration==60&(TE.Miss==1|TE.CorrectRejection==1)) + TE.StimulusOn(TE.StimulusDuration==60&(TE.Miss==1|TE.CorrectRejection==1));
stimon = TE.TrialStart(TE.StimulusDuration==60&(TE.Hit==1|TE.FalseAlarm==1)) + TE.StimulusOn(TE.StimulusDuration==60&(TE.Hit==1|TE.FalseAlarm==1));
stimon = TE.TrialStart(TE.StimulusDuration==60&TE.Hit==1) + TE.StimulusOn(TE.StimulusDuration==60&TE.Hit==1);
stimon = TE.TrialStart(TE.StimulusDuration==60) + TE.StimulusOn(TE.StimulusDuration==60);
stimon = TE.TrialStart(TE.StimulusDuration==60) + TE.LeftPortIn(TE.StimulusDuration==60);
stimon = TE.TrialStart(TE.FalseAlarm==1) + TE.StimulusOn(TE.FalseAlarm==1);
stimon = TE.TrialStart(TE.Miss==1) + TE.StimulusOn(TE.Miss==1);
stimon = TE.TrialStart(TE.CorrectRejection==1) + TE.StimulusOn(TE.CorrectRejection==1);
stimon = TE.TrialStart(TE.Hit==1) + TE.DeliverFeedback(TE.Hit==1);
stimon = TE.TrialStart(TE.Hit==1) + TE.LeftWaterValveOn(TE.Hit==1);
stimon = TE.TrialStart(TE.FalseAlarm==1) + TE.LeftPortIn(TE.FalseAlarm==1);
stimon = SE.BurstOn(~isnan(SE.BurstOn));
vdisc = round((stimon-starttime)*sr) + 1;

% STA
[st sta stats H] = stacall(vdisc,lfp,sr,wn);   % call STA code