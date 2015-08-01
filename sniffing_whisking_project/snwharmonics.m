function snwharmonics
%SNWHARMONICS   Coupling mode statistics.
%   SNWHARMONICS performs bootstrap test (1000 repetitions) to compare
%   sniffing frequency to the second harmonics of whisking frequency in the
%   1:2 mode and whisking frequency to the second harmonics of sniffing
%   frequency in the 2:1 mode. Period length distributions are calculated
%   from the intervals between inhalation start points (sniffing) or
%   detected whisking bouts (whisking). Frequency is estimated as the
%   reciprocal of the mean period length.
%
%   See also SNWMODESWITCH.

% Add Chronux to path
chronux_dir = fullfile(matlabroot,'work\chronux');
addpath(genpath(chronux_dir))

% Sampling rate
sr = 1000;
plotting = false;

% Load all segments
global DATADIR
fn = [DATADIR 'HSW\allsegments_speed_filtered.mat'];
load(fn)   % load all segments
fn = [DATADIR 'HSW\switches.mat'];
load(fn)   % load segment groups

% Progress indicator
wb = waitbar(0,'Please wait...','Name','Running SNWHARMONICS...');  % progress indicator
global WB
WB(end+1) = wb;

% Get raw data and calculate spectral variables
NumSeg = length(gr05);  % 1:2 mode
IFsniff = [];
IFwhisk = [];
for iS = 1:NumSeg   % loop through 1:2 segments
    cseg = allsegments(gr05(iS));   % current segment
    
    % Instanteous period length
    IFsniff = [IFsniff diff(cseg.inhalation_start)];   % sniffing period length
    IFwhisk = [IFwhisk diff(cseg.whisk_events)];   % whisking period length
    
    waitbar(iS/NumSeg)   % update progress indicator
end
close(wb)   % close progress indicator

% Mean frequencies
mn_sniff = 1 / mean(IFsniff);   % mean freq. for sniffing
mn_whisk = 1 / mean(IFwhisk);   % mean freq. for whisking

% Bootstrap confidence interval
bst = 1000;   % size of bootstrap sample
bootstrap_mn = nan(1,bst);
for k = 1:bst
    inx = randi(length(IFsniff),size(IFsniff));   % resampling
    bootstrap_mn(k) = 1 / mean(IFsniff(inx));   % bootstrap mean distribution
end
conf_int = prctile(bootstrap_mn,[2.5 97.5]);   % 5% confidence interval
p = sum(abs(2*mn_whisk-mn_sniff)<abs(bootstrap_mn-mn_sniff)) / bst;   % bootstrap test for difference between sniff freq. and second harmonic of whisk freq.

keyboard

% Progress indicator
wb = waitbar(0,'Please wait...','Name','Running SNWHARMONICS...');  % progress indicator
WB(end+1) = wb;

% Get raw data and calculate spectral variables
NumSeg = length(gr2);  % 2:1 mode
IFsniff = [];
IFwhisk = [];
for iS = 1:NumSeg   % loop through 2:1 segments
    cseg = allsegments(gr2(iS));   % current segment
    
    % Instanteous period length
    IFsniff = [IFsniff diff(cseg.inhalation_start)];   % sniffing period length
    IFwhisk = [IFwhisk diff(cseg.whisk_events)];   % whisking period length
    
    waitbar(iS/NumSeg)   % update progress indicator
end
close(wb)   % close progress indicator

% Mean frequencies
mn_sniff = 1 / mean(IFsniff);   % mean freq. for sniffing
mn_whisk = 1 / mean(IFwhisk);   % mean freq. for whisking

% Bootstrap confidence interval
bst = 1000;   % size of bootstrap sample
bootstrap_mn = nan(1,bst);
for k = 1:bst
    inx = randi(length(IFwhisk),size(IFwhisk));   % resampling
    bootstrap_mn(k) = 1 / mean(IFwhisk(inx));   % bootstrap mean distribution
end
conf_int = prctile(bootstrap_mn,[2.5 97.5]);   % 5% confidence interval
p = sum(abs(2*mn_sniff-mn_whisk)<abs(bootstrap_mn-mn_whisk)) / bst;   % bootstrap test for difference between sniff freq. and second harmonic of whisk freq.

keyboard