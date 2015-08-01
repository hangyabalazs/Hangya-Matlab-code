function somccgclust2
%SOMCCGCLUST2   Cluster analysis on response profiles.
%   SOMCCGCLUST2 performs hierarchical cluster analysis on PSTHs aligned to
%   3 events: HomeZoneIn, HomeZoneOut and WaterValveOn (output of
%   SOMRESPONSEPROFILES). Response profiles (PSTHs) are standardized and
%   clustered by GAP_STATISTICS, which tests whether an optimal number of
%   clusters results in a natural clustering of the data (see details
%   therein).
%
%   See also SOMRESPONSEPROFILES and GAP_STATISTICS.

% Directories
global DATADIR
global DATAPATH
inpdir = [DATAPATH 'SOM\ResponseCCG3\'];

% Import
load([DATADIR 'SOM_Sachin\PSTH\tagged_list_agreed'],'som_beh','pv_beh');
loadcb
NumCells = length(CELLIDLIST);
hzin = nan(NumCells,20);
hzout = nan(NumCells,20);
wvon = nan(NumCells,20);
for iC = 1:NumCells
    cellid = CELLIDLIST{iC};
    try
        fnm = [inpdir regexprep(cellid,'\.','_') '_HomeZoneIn1_CCG.mat'];
        ccs = load(fnm);
        hzin(iC,:) = standardize(ccs.ccr);   % z-score
        fnm = [inpdir regexprep(cellid,'\.','_') '_HomeZoneOut1_CCG.mat'];
        ccs = load(fnm);
        hzout(iC,:) = standardize(ccs.ccr);   % z-score
        fnm = [inpdir regexprep(cellid,'\.','_') '_WaterValveOn_CCG.mat'];
        ccs = load(fnm);
        wvon(iC,:) = standardize(ccs.ccr);   % z-score
    catch ME    % TrialEvents file could be missing
        if ~strncmp('Unable to read file',ME.message,19)
            error(ME.message)
        end
    end
end

% Remove NaNs
cellist = 1:NumCells;
inx = [];
for k = 1:NumCells
    if all(isnan(hzin(k,:))) || all(isnan(hzout(k,:))) || all(isnan(wvon(k,:)))
        inx = [inx k];
    end
end
hzin(inx,:) = [];
hzout(inx,:) = [];
wvon(inx,:) = [];
cellist(inx) = [];
cellids = CELLIDLIST;
cellids(inx) = [];

% Clustering
dec = 100;   % number of clusters
% [jnk gr_pv] = intersect(cellids,pv_beh);
% [jnk gr_som] = intersect(cellids,som_beh);
ccgS = nan([size(hzin) 3]);
ccgS(:,:,1) = hzin;
ccgS(:,:,2) = hzout;
ccgS(:,:,3) = wvon;
khat = gap_statistics(ccgS,dec);

% Visualization
% figure
% for k = 1:dec
%     eval(['subplot(' num2str(dec) ',1,' num2str(k) ')'])
%     imagesc(ccgS(c==k,:))
%     set(gca,'XTick',[])
%     colorbar
% end
% 
% arrayfun(@(k)length(intersect(CELLIDLIST(cellist(c==k)),pv_beh)),(1:dec))
% arrayfun(@(k)length(intersect(CELLIDLIST(cellist(c==k)),som_beh)),(1:dec))
% arrayfun(@(k)length(cellist(c==k)),(1:dec))
% 
% 
% clustered_cellids = cell(1,dec);
% for k = 1:dec
%     clustered_cellids{k} = CELLIDLIST(cellist(c==k));
% end