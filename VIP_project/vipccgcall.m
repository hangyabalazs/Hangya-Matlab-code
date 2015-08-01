function vipccgcall(I)
%VIPCCGCALL   Calls CCG for non-tetrode pairs.
%   VIPCCGCALL calculates cross-correlations for all non-tetrode pairs in
%   cell-base. For details, see SOMCCG_CONF_FILTER.
%
%   VIPCCGCALL(I) calls CCG for pairs of cells defined by the index set I.
%
%   See also SOMCCG_CONF_FILTER and XCORR_WRAND_FILTER.

% Pass the control to the user in case of error
dbstop if error

% Directories
global DATAPATH
fs = filesep;
resdir = [DATAPATH 'VIP' fs 'CCG_A12' fs];

% Determine time window
sr = 1000;      % sampling rate
wn = 30 * sr / 1000;    % 2 * 30 ms window

% Load CellBase
loadcb

% Input argument check
numCells = length(CELLIDLIST);
if nargin < 1
    I = 1:numCells;
end

% Cell pairs
PairOfCells = cell(0,2);
for iC = I
    cellid = CELLIDLIST{iC};
    [nm ps] = nontetrodepairs(cellid);
    for k = 1:nm
        [ratname session tetrode1 unit1] = cellid2tags(cellid);   % prevent duplicates
        [ratname session tetrode2 unit2] = cellid2tags(ps(k));
        if tetrode2 > tetrode1
            PairOfCells(end+1,1:2) = {cellid ps{k}};
        end
    end
end
numPairs = size(PairOfCells,1);

% Cell pair loop for CCG
limit_spikes = [100 50000];   % include max 50000 spikes; calculate only if min 100 spikes
[CCR LCCR UCCR] = deal(zeros(numPairs,2*wn+1));
[MeanH0 SDH0] = deal(nan(numPairs,1));
for iP = 1:numPairs
    cell1 = PairOfCells{iP,1};
    cell2 = PairOfCells{iP,2};
    [tseg ncc1 ncc2] = deal([]);
%     try
%         tseg = findSegs3(cell1,'segfilter','prestim2','prepulseinterval',1);   % for A1 
% %         tseg = findSegs3(cell1);   % for mPFC
%         ncc1 = extractSegSpikes(cell1,tseg);
%         ncc2 = extractSegSpikes(cell2,tseg);
%     catch
%         disp('Could not extract the segment.');
%     end
    ncc1 = loadcb(cell1,'SPIKES');
    ncc2 = loadcb(cell2,'SPIKES');
    ncc1 = ncc1(:);   % convert to column vectors
    ncc2 = ncc2(:);
        
    if length(ncc1) > limit_spikes(2);      % crop if too long to avoid out of memory
        ncc1 = ncc1(1:limit_spikes(2));
    end
    if length(ncc2) > limit_spikes(2);
        ncc2 = ncc2(1:limit_spikes(2));
    end
    
    if length(ncc1) > limit_spikes(1) && length(ncc2) > limit_spikes(1)     % minimum 100 spikes
        [H1 ccr lwr upr rccg] = somccg_conf_filter(ncc1,ncc2,wn);    % 1->2
%         ncl1 = regexprep(cell1,'.','_');
%         ncl2 = regexprep(cell2,'.','_');
        ncl1 = cell1;
        ncl1(ncl1=='.') = '_';
        ncl2 = cell2;
        ncl2(ncl2=='.') = '_';
        fnm = ['CCG_' ncl1 '_' ncl2 '.fig'];
        saveas(H1,fullfile(resdir,fnm))
        close(H1)
        CCR(iP,:) = ccr;
        LCCR(iP,:) = lwr;
        UCCR(iP,:) = upr;
        MeanH0(iP) = mean(mean(rccg,2),1);
        SDH0(iP) = mean(std(rccg,[],2),1);
        disp('running......')
    end
end

% Save
fnm = 'CCG_matrices.mat';
save(fullfile(resdir,fnm),'CCR','LCCR','UCCR','MeanH0','SDH0','PairOfCells')
keyboard