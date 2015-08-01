function [DM LrC Pref Pref2 valid_channels] = lightcluster(cellid,varargin)
%LIGHTCLUSTER   Clustering of light-evoked spikes.
%   [D L] = LIGHTCLUSTER(CELLID) analyzes whether the light-evoked spikes
%   of a putative tagged cell (CELLID) are separated from other
%   light-evoked spikes. It returns cluster quality indices (Isolation
%   Distance, D; L-ratio, L) for the distance of light-evoked spikes of
%   CELLID compared with all light-evoked spikes (i.e. cluster quality
%   indices restricted to light-evoked spikes).
%
%   [D L PREF PREF2] = LIGHTCLUSTER(CELLID) returns preferential occurrance
%   of light-evoked spikes in the tagged cluster with respect to all spikes
%   (PREF) or with respect to out-of-cluster spikes (PREF2).
%
%   [D L PREF PREF2 VALID_CHANNELS] = LIGHTCLUSTER(CELLID) also returns
%   channel validity.
%
%   See also LRATIO.

%   Edit log: BH 5/9/12, 6/20/12 

% Input arguments
prs = inputParser;
addRequired(prs,'cellid',@iscellid)
addParamValue(prs,'stim_period',[],@isnumeric)   % start and end point for the interval of light-evoked spikes
addParamValue(prs,'feature_names',{'Energy'},@(s)iscell(s)|ischar(s))  % names of MClust features to use
parse(prs,cellid,varargin{:})
g = prs.Results;
if ischar(g.feature_names)
    g.feature_names = {g.feature_names};
end

% Parse cellID
[r s t] = cellid2tags(cellid);

% Load spikes from Ntt file.
Nttfn = cellid2fnames(cellid,'Ntt');
all_spikes = LoadTT_NeuralynxNT(Nttfn);
TIMEFACTOR = getpref('cellbase','timefactor');    % scaling factor to convert spike times into seconds
all_spikes = all_spikes * TIMEFACTOR;
spk = loadcb(cellid,'Spikes');

% Latency of stimulated spikes
if isempty(g.stim_period)
    [lim1 lim2] = findStimPeriod(cellid);   % find putative stimulated period
    if isnan(lim1) || isnan(lim2)
        lim1 = 0.001;   % no activation detected
        lim2 = 0.006;
    end
else
    lim1 = g.stim_period(1);
    lim2 = g.stim_period(2);
end

% Evoked spikes
tsegs_evoked = rel2abstimes(cellid,[lim1 lim2],'stim','PulseOn');   % convert period to epochs relative to pulses
selts_evoked = extractSegSpikes(all_spikes,tsegs_evoked);   % find putative stimualated spikes
% selts_evoked = unique(selts_evoked);  % in very few cases, a spike time is registered twice
[junk junk2 evoked_inx] = intersect(selts_evoked,all_spikes); %#ok<*ASGLU> % get indices for light-evoked spikes
if ~isequal(junk,selts_evoked)   % internal check for spike times
    error('lightcluster:SpikeTimeMismatch','Mismatch between saved spike times and Ntt time stamps.')
end
[spk_evoked inx] = intersect(selts_evoked,spk); % light-evoked spikes of the tagged cell
cinx = setdiff(1:length(selts_evoked),inx);   % complementing index set

% Feature matrix
valid_channels = check_channel_validity(cellid);   % valid channels
X = [];
for k = 1:length(g.feature_names)
    basename = [getpref('cellbase','cell_pattern') num2str(t)];
    propfn = [basename '_' g.feature_names{k}];   % name of feature file (e.g. TT1_Amplitude)
    sessionpath = cellid2fnames(cellid,'sess');
    propfn_path = [sessionpath filesep 'FD'];   % where the feature file can be found
    if ~isdir(propfn_path)
        propfn_path = sessionpath;
    end
    propfn_full = [propfn_path filesep propfn];   % full path of feature file
    try
        wf_prop = load([propfn_full '.fd'],'-mat');     % load feature file
    catch    %#ok<CTCH>
        disp('Calculating missing feature data.')       % calculate feature file if it was not found
        calculate_features(sessionpath,propfn_path,g.feature_names(k),basename,valid_channels)
        wf_prop = load([propfn_full '.fd'],'-mat');     % load feature file
    end
    wf_prop = wf_prop.FeatureData(evoked_inx,:);    % only light-evoked spikes
    
    if ~isequal(size(wf_prop,2),sum(valid_channels))
        wf_prop  = wf_prop(:,valid_channels);   % estimated and original valid_channels don't match
    end
    X = [X; wf_prop']; %#ok<AGROW>
end

% Isolation Distance
XC = X(:,inx);
D = mahal(X',XC');
DC = D(cinx);
sDC = sort(DC);
linx = length(inx);
if linx <= length(sDC)
    DM = sDC(length(inx));
else
    DM = Inf;   % more spikes in the cluster than outside of it
end

% L-ratio
if ~isempty(DC)
    df = size(X,1);
    LC = sum(1-chi2cdf(DC,df));
    nc = size(XC,2);
    LrC = LC / nc;
else
    LrC = Inf;  % no spikes outside of the cluster: multiunit
end

% Preference of light-spikes
clus_allspikes = length(spk) / length(all_spikes);  % ratio of no. of spikes in the tagged cluster to all spikes
cluslight_alllight = length(spk_evoked) / length(selts_evoked);  % ratio of no. of light-evoked spikes in the tagged cluster to all light_evoked spikes
Pref = cluslight_alllight / clus_allspikes;  % preferential occurrance of light-evoked spikes in the tagged cluster
clus_allspikes2 = length(spk) / (length(all_spikes) - length(spk));  % ratio of no. of spikes in and out of the tagged cluster
cluslight_alllight2 = length(spk_evoked) / (length(selts_evoked) - length(spk_evoked));  % ratio of no. of light-evoked spikes in and out of the tagged cluster
Pref2 = cluslight_alllight2 / clus_allspikes2;  % preferential occurrance of light-evoked spikes in the tagged cluster

% -------------------------------------------------------------------------
function calculate_features(sessionpath,propfn_path,feature_names,basename,valid_channels)

% Create MClust variables
global MClust_FDdn
global MClust_ChannelValidity
global MClust_NeuralLoadingFunction
global MClust_TText
global MClust_FDext
global MClust_TTdn
global MClust_TTfn

MClust_FDext = '.fd';
MClust_TText = '.ntt';
MClust_TTfn = basename;
[t1, t2] = strtok(fliplr(which('MClust')),filesep);
MClust_Directory = fliplr(t2);
MClust_FDdn = propfn_path;
MClust_TTdn = sessionpath;
MClust_ChannelValidity = valid_channels;
MClust_NeuralLoadingFunction = char([MClust_Directory 'LoadingEngines\LoadTT_NeuralynxNT']);

% Calculate features
CalculateFeatures(basename,feature_names)