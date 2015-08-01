function sompvccgcall2(cellids)
%SOMPVCCGCALL2   Calls CCG for non-tetrode pairs.
%   SOMPVCCGCALL2 calculates cross-correlations for all non-tetrode pairs in
%   cell-base. For details, see SOMCCG_CONF_FILTER.
%
%   SOMPVCCGCALL2(I) calls CCG for pairs of cells defined by the index set I.
%   SOMPVCCGCALL2(CELLIDS) calls CCG for pairs of cells defined by the list
%   of CellIDs (CELLIDS).
%
%   See also SOMCCG_CONF_FILTER and XCORR_WRAND_FILTER.

% Pass the control to the user in case of error
dbstop if error

% Directories
global DATAPATH
fs = filesep;
resdir = [DATAPATH 'SOM' fs 'nt_pairs' fs];  % results directory
fnm = 'CCG_matrices3.mat';   % filename for saving the result matrices

% Include only long enough segments
longsegments = true;
seglim = 300;

% Determine time window
sr = 1000;      % sampling rate
wn = 30 * sr / 1000;    % 2 * 30 ms window

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

% Cell pairs
PairOfCells = cell(0,2);
numCells = length(cellids);  % number of cells
for iC = 1:numCells
    cellid = cellids{iC};
    [nm ps] = nontetrodepairs(cellid);
    ps = intersect(cellids,ps);
    nm = length(ps);
    for k = 1:nm
        [ratname session tetrode1 unit1] = cellid2tags(cellid);   % prevent duplicates
        [ratname session tetrode2 unit2] = cellid2tags(ps(k));
        if tetrode2 > tetrode1
            PairOfCells(end+1,1:2) = {cellid ps{k}}; %#ok<AGROW>
        end
    end
end
numPairs = size(PairOfCells,1);

% Cell pair loop for CCG
wb = waitbar(0,'Please wait...','Name','Running SOMPVCCGCALL2...');  % progress indicator
global WB
WB(end+1) = wb;
limit_spikes = [100 50000];   % include max 50000 spikes; calculate only if min 100 spikes
CCR = zeros(numPairs,2*wn+1);
LCCR = zeros(numPairs,2*wn+1);
UCCR = zeros(numPairs,2*wn+1);
MeanH0 = nan(numPairs,1);
SDH0 = nan(numPairs,1);
SegmentLength = nan(numPairs,1);
for iP = 1:numPairs
    cell1 = PairOfCells{iP,1};
    cell2 = PairOfCells{iP,2};
    try
        tseg = findSegs2_temp01(cell1);  % find time segments
        ltseg = tseg(2,:) - tseg(1,:);  % length of the segments
        if longsegments
            seginx = find(ltseg==max(ltseg));
            tseg = tseg(:,seginx(1));
            ltseg = ltseg(seginx(1));
            if tseg < seglim
                continue
            end
        end
        SegmentLength(iP) = sum(ltseg);
        [ncc1,seltsind,selisi] = extractSegSpikes(cell1,tseg);   % find spikes in the time segments
        [ncc2,seltsind,selisi] = extractSegSpikes(cell2,tseg);
    catch
        disp('Could not extract the segment.');
    end
%     ncc1 = loadcb(cell1,'SPIKES');   % use all spikes
%     ncc2 = loadcb(cell2,'SPIKES');
        
    if length(ncc1) > limit_spikes(2);      % crop if too long to avoid out of memory
        ncc1 = ncc1(1:limit_spikes(2));
    end
    if length(ncc2) > limit_spikes(2);
        ncc2 = ncc2(1:limit_spikes(2));
    end
    
    if length(ncc1) > limit_spikes(1) && length(ncc2) > limit_spikes(1)     % minimum 100 spikes
        [H1 ccr lwr upr rccg] = somccg_conf_filter(ncc1,ncc2,wn);    % 1->2
%         ncl1 = regexprep(cell1,'\.','_');
%         ncl2 = regexprep(cell2,'\.','_');
%         fnm = ['CCG_' ncl1 '_' ncl2 '.fig'];
%         saveas(H1,fullfile(resdir,fnm))   % save CCG plot
        close(H1)
        CCR(iP,:) = ccr;   % cross-correlogram
        LCCR(iP,:) = lwr;  % lower significance limit
        UCCR(iP,:) = upr;  % upper significance limit
        MeanH0(iP) = mean(mean(rccg,2),1);   % surrogate mean
        SDH0(iP) = mean(std(rccg,[],2),1);   % surrogate SD
        disp(['Pair #' num2str(iP) ' / ' num2str(numPairs) ' done......'])
    end
    
    % Save
    if isequal(mod(iP,50),0)   % save after every 50 pairs to prevent data loss
        save(fullfile(resdir,fnm),'PairOfCells','CCR','LCCR','UCCR','MeanH0','SDH0','SegmentLength')
        disp('Autosave done.')
    end
    
    waitbar(iP/numPairs)
end
close(wb)   % eliminate progress indicator

% Save
save(fullfile(resdir,fnm),'PairOfCells','CCR','LCCR','UCCR','MeanH0','SDH0','SegmentLength')