
% Load
global DATADIR
global DATAPATH
load([DATADIR 'SOM_Sachin\PSTH\tagged_list_agreed'],'som_beh','pv_beh');
loadcb
load([DATAPATH 'SOM\ResponseCCG3\CCG_STRENGTH.mat'])

% Subsets
[is som_inx inx2] = intersect(CELLIDLIST,som_beh);
[is pv_inx inx2] = intersect(CELLIDLIST,pv_beh);
nt_inx = setdiff((1:length(CELLIDLIST)),[som_inx pv_inx]);

% Final strength measure
UaS = squeeze(CcgStrength(:,:,2));
LaS = squeeze(CcgStrength(:,:,3));
UaStrength = squeeze(CcgStrength(:,:,4));
for k = 1:size(CcgStrength,1)
    tuas = UaStrength(k,:);
%     tu = UaS(k,:);
%     tl = LaS(k,:);
%     tinx = tu<0.0001&tl<0.0001;
%     tuas(tinx) = NaN;
    UaStrength(k,:) = tuas / nansum(tuas);
end

% Remove NaNs
ccgS = UaStrength;
NumCells = size(CcgStrength,1);
cellist = 1:NumCells;
inx = [];
for k = 1:NumCells
    if all(isnan(ccgS(k,:))) || all(ccgS(k,:)==0)
        inx = [inx k];
    end
end
ccgS(inx,:) = [];
cellist(inx) = [];
cellids = CELLIDLIST;
cellids(inx) = [];

% Clustering
dec = 30;   % number of clusters
[jnk gr_pv] = intersect(cellids,pv_beh);
[jnk gr_som] = intersect(cellids,som_beh);
[c cc] = somclust(nan2zero(ccgS),dec,gr_pv,gr_som);

% Visualization
figure
for k = 1:dec
    eval(['subplot(' num2str(dec) ',1,' num2str(k) ')'])
    imagesc(ccgS(c==k,:))
    set(gca,'XTick',[])
    colorbar
end

arrayfun(@(k)length(intersect(CELLIDLIST(cellist(c==k)),pv_beh)),(1:dec))
arrayfun(@(k)length(intersect(CELLIDLIST(cellist(c==k)),som_beh)),(1:dec))
arrayfun(@(k)length(cellist(c==k)),(1:dec))


clustered_cellids = cell(1,dec);
for k = 1:dec
    clustered_cellids{k} = CELLIDLIST(cellist(c==k));
end