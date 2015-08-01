function spikeshapeanalysis(cellids,issave)
%SPIKESHAPEANALYSIS   Spike shape features.
%   SPIKESHPEANALYSIS(CELLIDS,ISSAVE) calculates a battery of spike
%   waveform features including magnitudes (e.g. peak, pre-valley,
%   post-valley), times (e.g. spike width = HHMW, peak-to-valley and
%   peak-to-baseline times in microseconds) and area features (e.g. full
%   integral, AHP). The results are stored in a cell array of structures.
%   Input parameters: 
%       CELLIDS - list of cell IDs or index set to CELLIDLIST (see CellBase
%           documentation); if empty or not specified, all cells are
%           selected from the CellBase
%       ISSAVE - controls saving
%
%   See also PLOTWAVEFORMS and EXTRACTSPIKEWAVEFORMS.

%   Balazs Hangya, Cold Spring Harbor Laboratory
%   1 Bungtown Road, Cold Spring Harbor
%   balazs.cshl@gmail.com
%   21-Oct-2013

%   Edit log: BH 21/10/13

% Pass the control to the user in case of error
dbstop if error

% Input argument check
error(nargchk(0,2,nargin))
if nargin < 2
    issave = true;
end
if nargin < 1 || isempty(cellids)
    loadcb   % load CellBase
    cellids = CELLIDLIST;
else
    if isnumeric(cellids)
        loadcb   % load CellBase
        cellids = CELLIDLIST(cellids);   % index set to CELLIDLIST
    elseif ischar(cellids)
        cellids = {cellids};   % only one cellID passed
    elseif iscellstr(cellids)
        % list of cell IDs
    else
        error('spikeshapeanalysis:inputArg','Unsupported format for cell IDs.')
    end
end
NumCells = length(cellids);   % number of cells

% Directories
global DATAPATH
resdir = fullfile(DATAPATH,'HDB','spikeshapeanalysis_newdata');

% Code kernel
SpikeShape = cell(1,NumCells);
fnmt = fullfile(resdir,'temp.mat');   % temporary file name
for k = 1:NumCells   % loop through cell IDs
    cellid = cellids{k};
    disp([num2str(k) '   ' cellid])
    SpikeShape{k} = main(cellid);   % calculate spike features
    
    if isequal(mod(k,20),0)   % save after every 20 cells to prevent data loss
        save(fnmt,'SpikeShape','cellids')
        disp('Autosave done.')
    end
end

% Save
if issave
    fnm = fullfile(resdir,'SpikeShape.mat');
    save(fnm,'SpikeShape','cellids')
end

% -------------------------------------------------------------------------
function SpikeShape = main(cellid)

% Average waveform from the largest channel
wv = extractSpikeWaveforms(cellid,'all','chans','mean_max');

% Waveform variables
sr = 32552;     % DigiLynx sampling rate
lwv = length(wv);   % number of data points
[mx maxtime] = max(wv);   % maximum and peak location
[mn mintime] = min(wv);   % minimum and trough location
[prevalley prevalleytime] = min(wv(1:maxtime));   % pre-valley
isprevalley = prevalley < 0;
[postvalley postvalleytime] = min(wv(maxtime+1:end));   % post-valley
ispostvalley = postvalley < 0;
postvalleytime = maxtime + postvalleytime;
zerotime = valuecrossing(1:lwv,wv,0);   % baseline crossing
prepeak0time = zerotime(find(zerotime<maxtime,1,'last'));   % baseline before peak
if isempty(prepeak0time)
    prepeak0time = 1;
end
postpeak0time =  zerotime(find(zerotime>maxtime,1,'first'));   % baseline after peak
if isempty(postpeak0time)
    postpeak0time = lwv;
end
prevalley0time =  zerotime(find(zerotime<prevalleytime,1,'last'));   % baseline before pre-valley
if isempty(prevalley0time)
    prevalley0time = 1;
end
postvalley0time =  zerotime(find(zerotime>postvalleytime,1,'first'));   % baseline after post-valley
if isempty(postvalley0time)
    postvalley0time = lwv;
end
HHtime = valuecrossing(1:lwv,wv,mx/2);   % half-hight crossing
HH1time = HHtime(find(HHtime<maxtime,1,'last'));   % guard against multiple crossings
HH2time = HHtime(find(HHtime>maxtime,1,'first'));
spikewidth = HH2time - HH1time;   % spike width (MWHH) in data points

% Waveform features - magnitude
SpikeShape.Spike = wv;   % save the waveform
SpikeShape.Peak = mx;   % peak magnitude
SpikeShape.Valley = abs(mn);   % valley magnitude
SpikeShape.Amplitude = mx - mn;   % full amplitude
SpikeShape.IsPreValley = isprevalley;   % if prevalley < 0
SpikeShape.PreValley = abs(prevalley);   % pre-valley magnitude
SpikeShape.IsPostValley = ispostvalley;   % if postvalley < 0
SpikeShape.PostValley = abs(postvalley);   % pre-valley magnitude
SpikeShape.PeakToPostValley = SpikeShape.Peak / SpikeShape.PostValley;   % peak-to-valley ratio
SpikeShape.PeakToPreValley = SpikeShape.Peak / SpikeShape.PreValley;   % peak-to-prevalley ratio
SpikeShape.PostValleyToPeak = SpikeShape.PostValley / SpikeShape.Peak;   % post-valley normalized to peak
SpikeShape.PreValleyToPeak = SpikeShape.PreValley / SpikeShape.Peak;   % pre-valley normalized to peak

% Waveform features - time
SpikeShape.SpikeWidth = spikewidth / sr * 1e6;   % spike width in microseconds
SpikeShape.PeakTime = (maxtime - 1) / sr * 1e6;   % peak time
SpikeShape.PrePeak0Time = (prepeak0time - 1) / sr * 1e6;   % pre-peak time
SpikeShape.PostPeak0Time = (postpeak0time - 1) / sr * 1e6;   % post-peak time
SpikeShape.PreValleyTime = (prevalleytime - 1) / sr * 1e6;   % pre-valley time
SpikeShape.PostValleyTime = (postvalleytime - 1) / sr * 1e6;   % post-valley time
SpikeShape.PreValley0Time = (prevalley0time - 1) / sr * 1e6;   % pre-valley time
SpikeShape.PostValley0Time = (postvalley0time - 1) / sr * 1e6;   % post-valley time
SpikeShape.PreValleyToPeakTime = (maxtime - prevalleytime) / sr * 1e6;  % pre-valley to peak time in microseconds
SpikeShape.PeakToPostValleyTime = (postvalleytime - maxtime) / sr * 1e6;  % peak to post-valley time in microseconds
SpikeShape.BaselineToPeakTime = (maxtime - prepeak0time) / sr * 1e6;  % baseline to peak time in microseconds
SpikeShape.PeakToBaselineTime = (postpeak0time - maxtime) / sr * 1e6;  % peak to baseline time in microseconds
SpikeShape.PrePeakTime = (maxtime - prevalley0time) / sr * 1e6;  % first baseline to peak time in microseconds
SpikeShape.PostPeakTime = (postvalley0time - maxtime) / sr * 1e6;  % peak to last baseline time in microseconds

% Waveform features - area
SpikeShape.Integral = sum(abs(wv));   % spike integral
SpikeShape.AHP = sum(abs(wv(ceil(postpeak0time):floor(postvalley0time))));   % after-hyperpolarization
SpikeShape.BHP = sum(abs(wv(ceil(prevalley0time):floor(prepeak0time))));   % before-hyperpolarization
SpikeShape.PeakIntegral = sum(abs(wv(ceil(prepeak0time):floor(postpeak0time))));   % peak integral
SpikeShape.HalfPeakIntegral = sum(abs(wv(ceil(HH1time):floor(HH2time))));   % peak integral between HH crossings
SpikeShape.AHPToIntegral = SpikeShape.AHP / SpikeShape.Integral;   % AHP normalized to integral
SpikeShape.BHPToIntegral = SpikeShape.BHP / SpikeShape.Integral;   % BHP normalized to integral
SpikeShape.PeakIntegralToIntegral = SpikeShape.PeakIntegral / SpikeShape.Integral;   % peak integral normalized to integral
SpikeShape.HalfPeakIntegralToIntegral = SpikeShape.HalfPeakIntegral / SpikeShape.Integral;   % half peak integral normalized to integral