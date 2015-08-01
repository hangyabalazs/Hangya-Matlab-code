function ccg(cellids,varargin)
%CCG   Cross-correlation.
%   CCG calculates cross-correlations. Window size is set to +-1000 ms with
%   a 1 ms resolution. Maximum 50000 spikes are included to avoid memory
%   problems. CCG is not calculated if one of the cells has less than 100
%   spikes. Segments are filtered with FINDSEGS3. Minimal shift for
%   shuffled CCGs is set to 1100 ms. For details on the algorithm, see
%   SOMCCG_CONF_FILTER.
%
%   CCG(I) calls CCG for pairs of cells defined by the cell ID list 
%   (or index set to CELLIDLIST, see CellBase documentation) I. By default,
%   non-tetrode pairs within the given list are selected for analysis.
%   Optional input arguments (parameter-value pairs with default values):
%       'issave', false - controls saving behavior; plots and
%           cross-correlation matrices with confidence intervals are saved 
%           only if 'issave' is set to true
%       'whichcells', 'nontetrodepairs' - method of pair selection;
%           'nontetrodepairs' selects cells from other tetrodes,
%           'tetrodepairs' selects cells from the same tetrode and
%           'allpairs' selects all cells from the session
%       'include', 'list' - by default, only pairs for which both cells are
%           included in I are analyzed; if 'include' is set to 'cellbase',
%           all cells in CellBase that are paired with the ones in I
%           according to 'whichcells' are analyzed
%
%   See also ACG, SOMCCG_CONF_FILTER and XCORR_WRAND_FILTER.

%   Balazs Hangya, Cold Spring Harbor Laboratory
%   1 Bungtown Road, Cold Spring Harbor
%   balazs.cshl@gmail.com
%   14-June-2013

% Input arguments
prs = inputParser;
addRequired(prs,'cellids',@(s)iscell(s)|iscellstr(s)|ischar(s))
addParamValue(prs,'issave',false,@islogical)   % control saving behavior
addParamValue(prs,'whichcells','nontetrodepairs',...
    @(s)ismember(s,{'nonterodepairs','tetrodepairs','allpairs'}))   % which cells to include
addParamValue(prs,'include','list',...
    @(s)ismember(s,{'list','cellbase'}))   % cell pair selection behavior: only from input list or full CellBase
parse(prs,cellids,varargin{:})
g = prs.Results;

% Pass the control to the user in case of error
dbstop if error

% Directories
global DATAPATH
fs = filesep;
resdir = [DATAPATH 'NB\CCG2\_\'];  % results directory
fnmm = 'CCG_matrices.mat';   % filename for saving the result matrices

% Include only long enough segments
longsegments = false;  % control whether to use this option
seglim = 300;

% Determine time window
sr = 1000;      % sampling rate
wn = 1000 * sr / 1000;    % 2 * 1000 ms window

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
    [animalID sessionID tetrode1 unit1] = cellid2tags(cellid);
    switch g.whichcells
        case 'nontetrodepairs'
            [nm ps] = nontetrodepairs(cellid);   % cells on other tetrodes
        case 'tetrodepairs'
            [nm ps] = tetrodepairs(cellid);   % cells on the same tetrode
            ps = setdiff(ps,cellid);
            nm = length(ps);
        case 'allpairs'
            ps = findcell('rat',animalID,'session',sessionID);   % all concurrently recorded cells
            ps = setdiff(ps,cellid);
            nm = length(ps);
    end
    if isequal(g.include,'list')
        ps = intersect(cellids,ps);   % include only those that are in the input list
        nm = length(ps);
    end
    for k = 1:nm
        [j1, j2, tetrode2, unit2] = cellid2tags(ps(k));
        if (tetrode2*10+unit2) > (tetrode1*10+unit1)   % prevent duplicates
            PairOfCells(end+1,1:2) = {cellid ps{k}}; %#ok<AGROW>
        end
    end
end
numPairs = size(PairOfCells,1);

% Cell pair loop for CCG
wb = waitbar(0,'Please wait...','Name','Running CCG...');  % progress indicator
global WB
WB(end+1) = wb;
limit_spikes = [100 50000];   % include max 50000 spikes; calculate only if min 100 spikes
[CCR LCCR UCCR] = deal(zeros(numPairs,2*wn+1));
[MeanH0 SDH0] = deal(nan(numPairs,1));
SDH0 = nan(numPairs,1);
SegmentLength = nan(numPairs,1);
for iP = 1:numPairs   % loop through pairs of cells
    cell1 = PairOfCells{iP,1};
    cell2 = PairOfCells{iP,2};
    try
%         tseg = findSegs3(cell1,'segfilter','stimfb_excl_nb',...
%             'light_activation_duration',[-5 5],'margins',[0 0]);  % find time segments
%         tseg = findSegs3(cell1,'segfilter','prestim3');  % find time segments
        tseg = findSegs3(cell1,'segfilter','fb_incl_nb',...
            'feedback_duration',[1.5 2.5],'margins',[0 0],'min_int',0);  % find time segments
%         tseg = findSegs3(cell1,'segfilter','cue_incl_nb',...
%             'feedback_duration',[-1.5 0],'margins',[0 0],'min_int',0);  % find time segments
        ltseg = tseg(2,:) - tseg(1,:);  % length of the segments
        if longsegments   % later to be implemented as an input option
            seginx = find(ltseg==max(ltseg));
            tseg = tseg(:,seginx(1));   % find the longest segment
            ltseg = ltseg(seginx(1));
            if tseg < seglim   % use the longest segment if it's longer than the threshold
                continue
            end
        end
        SegmentLength(iP) = sum(ltseg);  % cumulative length of the segments
        [ncc1,seltsind,selisi] = extractSegSpikes(cell1,tseg);   % find spikes in the time segments
        [ncc2,seltsind,selisi] = extractSegSpikes(cell2,tseg);
    catch ME
        disp(ME.message)
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
        [H1 ccr lwr upr rccg] = somccg_conf_filter(ncc1,ncc2,wn,1100);    % 1->2
        xl = xlim;   % put the cell IDs in the plot
        yl = ylim;
        text(xl(1)+(xl(2)-xl(1))*0.7,yl(1)+(yl(2)-yl(1))*0.9,regexprep(cell1,'_',' '))
        text(xl(1)+(xl(2)-xl(1))*0.7,yl(1)+(yl(2)-yl(1))*0.85,regexprep(cell2,'_',' '))
        if g.issave   % save figure
            ncl1 = regexprep(cell1,'\.','_');
            ncl2 = regexprep(cell2,'\.','_');
            fnm = ['CCG_' ncl1 '_' ncl2 '.fig'];
            saveas(H1,fullfile(resdir,fnm))   % save CCG plot
            close(H1)
        end
        
        ccr2 = reshape(ccr(1:end-1),10,(length(ccr)-1)/10);
        sccr2 = sum(ccr2);
        nqf = sr / 10 / 2;
        flt = fir1(64,4/nqf,'high');
        fccr = filtfilt(flt,1,sccr2);
        figure
        bar(linspace(-wn,wn,length(fccr)),fccr,'BarWidth',1,...
            'EdgeColor','k','FaceColor','k')
        
        CCR(iP,:) = ccr;   % cross-correlogram
        LCCR(iP,:) = lwr;  % lower significance limit
        UCCR(iP,:) = upr;  % upper significance limit
        MeanH0(iP) = mean(mean(rccg,2),1);   % surrogate mean
        SDH0(iP) = mean(std(rccg,[],2),1);   % surrogate SD
        disp(['Pair #' num2str(iP) ' / ' num2str(numPairs) ' done......'])
    end
    
    % Save
    if g.issave
        if isequal(mod(iP,50),0)   % save after every 50 pairs to prevent data loss
            save(fullfile(resdir,fnm),'PairOfCells','CCR','LCCR','UCCR','MeanH0','SDH0','SegmentLength')
            disp('Autosave done.')
        end
    end
    
    waitbar(iP/numPairs)
end
close(wb)   % eliminate progress indicator

% Save
if g.issave
    save(fullfile(resdir,fnmm),'PairOfCells','CCR','LCCR','UCCR','MeanH0','SDH0','SegmentLength')
end