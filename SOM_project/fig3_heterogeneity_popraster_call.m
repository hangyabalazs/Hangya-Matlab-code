%% load

load('C:\Balazs\_data\SOM_Sachin\clustered_cellids\clustered_cellids.mat')
sel_clus_inx = [6 81-3 35-3 5 18 7 2 16 17 48 4 23 39-11 2 6 32 7 34-8 9 60 26 7 22-6 13 46 11 8 17 29-13];
sel_clus_inx = [sel_clus_inx(1:26) NaN sel_clus_inx(27:29)];
% shuff_inx = [7 4 11 15 20 22 10 23 18 17 9 24 1 26 25 19 27 6 3 21 28 8 16 5 29 2 14 13 12];

%% cellids

cellids = {};
for k = [1:26 28:30]
%     cellids(end+1) = clustered_cellids{k}(6);
    cellids(end+1) = clustered_cellids{k}(sel_clus_inx(k));
end
cellids(1) = clustered_cellids{1}(50);
cellids(2) = clustered_cellids{2}(17);
cellids(5) = clustered_cellids{5}(26);
cellids(9) = clustered_cellids{9}(23);
cellids(15) = clustered_cellids{15}(3);
cellids(24) = clustered_cellids{24}(20);
cellids(25) = clustered_cellids{25}(30);
cellids(26) = clustered_cellids{26}(50);
cellids(28) = clustered_cellids{29}(19);
cellids(27) = clustered_cellids{28}(16);
% [Y,I] = sort(shuff_inx);
% cellids = cellids(I);

%% cellids for a cluster
% cellids = clustered_cellids{9};

%% PSTH

R2 = fig3_heterogeneity_popraster(cellids);
[Y,I] = max(R2,[],2);
[Y,I] = sort(I);
R2 = fig3_heterogeneity_popraster(cellids(I));
