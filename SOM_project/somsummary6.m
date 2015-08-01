%% load

load('C:\My Dropbox\KepecsLab\_Duda\for balazs\pop_psth_for_PCA\HomeZout.mat')
pv_MatrixPsth1 = pv;
som_MatrixPsth1 = som;
non_tagged_MatrixPsth1 = non_tagged_WS;
Pv_neighbors_WS1 = Pv_neighbors_WS;
Som_neighbors_WS1 = Som_neighbors_WS;

load('C:\My Dropbox\KepecsLab\_Duda\for balazs\pop_psth_for_PCA\HomeZin.mat')
pv_MatrixPsth2 = pv;
som_MatrixPsth2 = som;
non_tagged_MatrixPsth2 = non_tagged_WS;
Pv_neighbors_WS2 = Pv_neighbors_WS;
Som_neighbors_WS2 = Som_neighbors_WS;

%% remove NaNs

[x1 y1] = find(isnan(non_tagged_MatrixPsth1));
[x2 y2] = find(isnan(non_tagged_MatrixPsth2));
x = union(unique(x1),unique(x2));
non_tagged_MatrixPsth1(x,:) = [];
non_tagged_MatrixPsth2(x,:) = [];

%% restrict in time - new data

pv_MatrixPsth1 = pv_MatrixPsth1(:,7:207);
som_MatrixPsth1 = som_MatrixPsth1(:,7:207);
non_tagged_MatrixPsth1 = non_tagged_MatrixPsth1(:,7:207);
Pv_neighbors_WS1 = Pv_neighbors_WS1(:,7:207);
Som_neighbors_WS1 = Som_neighbors_WS1(:,7:207);
pv_MatrixPsth2 = pv_MatrixPsth2(:,7:207);
som_MatrixPsth2 = som_MatrixPsth2(:,7:207);
non_tagged_MatrixPsth2 = non_tagged_MatrixPsth2(:,7:207);
Pv_neighbors_WS2 = Pv_neighbors_WS2(:,7:207);
Som_neighbors_WS2 = Som_neighbors_WS2(:,7:207);

%% smooth

sizepv = size(pv_MatrixPsth1);
spv1 = zeros(sizepv(1),sizepv(2));
for k = 1:sizepv(1)
    spv1(k,:) = smooth(pv_MatrixPsth1(k,:),'linear',11);
end

sizesom = size(som_MatrixPsth1);
ssom1 = zeros(sizesom(1),sizesom(2));
for k = 1:sizesom(1)
    ssom1(k,:) = smooth(som_MatrixPsth1(k,:),'linear',11);
end

sizent = size(non_tagged_MatrixPsth1);
snt1 = zeros(sizent(1),sizent(2));
for k = 1:sizent(1)
    snt1(k,:) = smooth(non_tagged_MatrixPsth1(k,:),'linear',11);
end

pvn1 = zeros(size(Pv_neighbors_WS1));
for k = 1:size(Pv_neighbors_WS1,1)
    pvn1(k,:) = smooth(Pv_neighbors_WS1(k,:),'linear',11);
end

somn1 = zeros(size(Som_neighbors_WS1));
for k = 1:size(Som_neighbors_WS1,1)
    somn1(k,:) = smooth(Som_neighbors_WS1(k,:),'linear',11);
end

spv2 = zeros(sizepv(1),sizepv(2));
for k = 1:sizepv(1)
    spv2(k,:) = smooth(pv_MatrixPsth2(k,:),'linear',11);
end

ssom2 = zeros(sizesom(1),sizesom(2));
for k = 1:sizesom(1)
    ssom2(k,:) = smooth(som_MatrixPsth2(k,:),'linear',11);
end

snt2 = zeros(sizent(1),sizent(2));
for k = 1:sizent(1)
    snt2(k,:) = smooth(non_tagged_MatrixPsth2(k,:),'linear',11);
end

pvn2 = zeros(size(Pv_neighbors_WS2));
for k = 1:size(Pv_neighbors_WS2,1)
    pvn2(k,:) = smooth(Pv_neighbors_WS2(k,:),'linear',11);
end

somn2 = zeros(size(Som_neighbors_WS2));
for k = 1:size(Som_neighbors_WS2,1)
    somn2(k,:) = smooth(Som_neighbors_WS2(k,:),'linear',11);
end

%% supression index

template = mean(spv1);
spvd2 = zeros(1,sizepv(1));
for k = 1:sizepv(1)
    spvd2(k) = sum(spv1(k,:).*template);
end

pvnd2 = zeros(1,size(pvn1,1));
for k = 1:size(pvn1,1)
    pvnd2(k) = sum(pvn1(k,:).*template);
end

ssomd2 = zeros(1,sizesom(1));
for k = 1:sizesom(1)
    ssomd2(k) = sum(ssom1(k,:).*template);
end

% inx1 = 125:150;
% inx2 = 1:100;
inx2 = 111:200;
inx1 = 80:110;
si_pv2 = zeros(1,sizepv(1));
for k = 1:sizepv(1)
    si_pv2(k) = mean(spv2(k,inx2)) - mean(spv2(k,inx1));
end

si_som2 = zeros(1,sizesom(1));
for k = 1:sizesom(1)
    si_som2(k) = mean(ssom2(k,inx2)) - mean(ssom2(k,inx1));
end

si_pvn2 = zeros(1,size(pvn2,1));
for k = 1:size(pvn2,1)
    si_pvn2(k) = mean(pvn2(k,inx2)) - mean(pvn2(k,inx1));
end

si_somn2 = zeros(1,size(somn2,1));
for k = 1:size(somn2,1)
    si_somn2(k) = mean(somn2(k,inx2)) - mean(somn2(k,inx1));
end

si_snt2 = zeros(1,size(snt2,1));
for k = 1:size(snt2,1)
    si_snt2(k) = mean(snt2(k,inx2)) - mean(snt2(k,inx1));
end

figure
plot(spvd2,si_pv2,'r.','MarkerSize',20)
hold on
plot(pvnd2,si_pvn2,'y.','MarkerSize',20)
plot(ssomd2,si_som2,'b.','MarkerSize',20)

figure
[nm xout] = hist(si_som2);
stairs(xout,nm/sum(nm),'Color','b','LineWidth',3)
hold on
[nm xout] = hist(si_pv2);
stairs(xout,nm/sum(nm),'Color','r')
[nm xout] = hist(si_snt2);
stairs(xout,nm/sum(nm),'Color',[0.7 0.7 0.7])

%% templates again

template = mean(spv1);
spvd1 = zeros(1,sizepv(1));
for k = 1:sizepv(1)
    spvd1(k) = sum(spv1(k,:).*template);
end

ssomd1 = zeros(1,sizesom(1));
for k = 1:sizesom(1)
    ssomd1(k) = sum(ssom1(k,:).*template);
end

pvnd1 = zeros(1,size(pvn1,1));
for k = 1:size(pvn1,1)
    pvnd1(k) = sum(pvn1(k,:).*template);
end

somnd1 = zeros(1,size(somn1,1));
for k = 1:size(somn1,1)
    somnd1(k) = sum(somn1(k,:).*template);
end

sntd1 = zeros(1,size(snt1,1));
for k = 1:size(snt1,1)
    sntd1(k) = sum(snt1(k,:).*template);
end

template2 = mean(ssom1);
spve1 = zeros(1,sizepv(1));
for k = 1:sizepv(1)
    spve1(k) = sum(spv1(k,:).*template2);
end

ssome1 = zeros(1,sizesom(1));
for k = 1:sizesom(1)
    ssome1(k) = sum(ssom1(k,:).*template2);
end

pvne1 = zeros(1,size(pvn1,1));
for k = 1:size(pvn1,1)
    pvne1(k) = sum(pvn1(k,:).*template2);
end

somne1 = zeros(1,size(somn1,1));
for k = 1:size(somn1,1)
    somne1(k) = sum(somn1(k,:).*template2);
end

snte1 = zeros(1,size(snt1,1));
for k = 1:size(snt1,1)
    snte1(k) = sum(snt1(k,:).*template2);
end

figure
plot(spvd1,spve1,'r.','MarkerSize',20)
hold on
plot(ssomd1,ssome1,'b.','MarkerSize',20)
plot(sntd1,snte1,'k.','MarkerSize',20)

template = mean(spv2);
spvd2 = zeros(1,sizepv(1));
for k = 1:sizepv(1)
    spvd2(k) = sum(spv2(k,:).*template);
end

ssomd2 = zeros(1,sizesom(1));
for k = 1:sizesom(1)
    ssomd2(k) = sum(ssom2(k,:).*template);
end

pvnd2 = zeros(1,size(pvn1,1));
for k = 1:size(pvn1,1)
    pvnd2(k) = sum(pvn2(k,:).*template);
end

somnd2 = zeros(1,size(somn1,1));
for k = 1:size(somn1,1)
    somnd2(k) = sum(somn2(k,:).*template);
end

sntd2 = zeros(1,size(snt2,1));
for k = 1:size(snt2,1)
    sntd2(k) = sum(snt2(k,:).*template);
end

template2 = mean(ssom2);
spve2 = zeros(1,sizepv(1));
for k = 1:sizepv(1)
    spve2(k) = sum(spv2(k,:).*template2);
end

ssome2 = zeros(1,sizesom(1));
for k = 1:sizesom(1)
    ssome2(k) = sum(ssom2(k,:).*template2);
end

pvne2 = zeros(1,size(pvn1,1));
for k = 1:size(pvn1,1)
    pvne2(k) = sum(pvn2(k,:).*template2);
end

somne2 = zeros(1,size(somn1,1));
for k = 1:size(somn1,1)
    somne2(k) = sum(somn2(k,:).*template2);
end

snte2 = zeros(1,size(snt2,1));
for k = 1:size(snt2,1)
    snte2(k) = sum(snt2(k,:).*template2);
end

figure
plot(spvd2,spve2,'r.','MarkerSize',20)
hold on
plot(ssomd2,ssome2,'b.','MarkerSize',20)
plot(sntd2,snte2,'k.','MarkerSize',20)

%% supression time

template = mean(spv1);
spvd = zeros(1,sizepv(1));
for k = 1:sizepv(1)
    spvd(k) = sum(spv1(k,:).*template);
end

pvnd = zeros(1,size(pvn1,1));
for k = 1:size(pvn1,1)
    pvnd(k) = sum(pvn1(k,:).*template);
end

ssomd = zeros(1,sizesom(1));
for k = 1:sizesom(1)
    ssomd(k) = sum(ssom1(k,:).*template);
end

somnd = zeros(1,size(somn1,1));
for k = 1:size(somn1,1)
    somnd(k) = sum(somn1(k,:).*template);
end

ntd = zeros(1,size(snt1,1));
for k = 1:size(snt1,1)
    ntd(k) = sum(snt1(k,:).*template);
end

mlt = 0.25;   % 0.25: maybe best; 0.3 also good
p1 = 1;
p2 = 99;
los_pv1 = zeros(1,size(spv1,1));
bs = (max(spv1(:)) - min(spv1(:))) * mlt + min(spv1(:));
bs = (prctile(spv1(:),p2) - prctile(spv1(:),p1)) * mlt + prctile(spv1(:),p1);
for k = 1:size(spv1,1)
    los_pv1(k) = sum(spv1(k,:)<bs);
end

los_som1 = zeros(1,size(ssom1,1));
bs = (max(ssom1(:)) - min(ssom1(:))) * mlt + min(ssom1(:));
bs = (prctile(ssom1(:),p2) - prctile(ssom1(:),p1)) * mlt + prctile(ssom1(:),p1);
for k = 1:size(ssom1,1)
    los_som1(k) = sum(ssom1(k,:)<bs);
end

los_pvn1 = zeros(1,size(pvn1,1));
bs = (max(pvn1(:)) - min(pvn1(:))) * mlt + min(pvn1(:));
bs = (prctile(pvn1(:),p2) - prctile(pvn1(:),p1)) * mlt + prctile(pvn1(:),p1);
for k = 1:size(pvn1,1)
    los_pvn1(k) = sum(pvn1(k,:)<bs);
end

los_somn1 = zeros(1,size(somn1,1));
bs = (max(somn1(:)) - min(somn1(:))) * mlt + min(somn1(:));
bs = (prctile(somn1(:),p2) - prctile(somn1(:),p1)) * mlt + prctile(somn1(:),p1);
for k = 1:size(somn1,1)
    los_somn1(k) = sum(somn1(k,:)<bs);
end

los_nt1 = zeros(1,size(snt1,1));
bs = (max(snt1(:)) - min(snt1(:))) * mlt + min(snt1(:));
bs = (prctile(snt1(:),p2) - prctile(snt1(:),p1)) * mlt + prctile(snt1(:),p1);
for k = 1:size(snt1,1)
    los_nt1(k) = sum(snt1(k,:)<bs);
end

figure
plot(ntd,los_nt1,'.','MarkerSize',20,'Color',[0.7 0.7 0.7])
hold on
plot(pvnd,los_pvn1,'y.','MarkerSize',20)
plot(spvd,los_pv1,'r.','MarkerSize',20)
plot(ssomd,los_som1,'b.','MarkerSize',20)

%% supression time #2

template = mean(ssom2);
spvd = zeros(1,sizepv(1));
for k = 1:sizepv(1)
    spvd(k) = sum(spv2(k,:).*template);
end

somnd = zeros(1,size(somn2,1));
for k = 1:size(somn2,1)
    somnd(k) = sum(somn2(k,:).*template);
end

ssomd = zeros(1,sizesom(1));
for k = 1:sizesom(1)
    ssomd(k) = sum(ssom2(k,:).*template);
end

pvnd = zeros(1,size(pvn2,1));
for k = 1:size(pvn2,1)
    pvnd(k) = sum(pvn2(k,:).*template);
end

mlt = 0.25;   % 0.25: maybe best; 0.3 also good
p1 = 1;
p2 = 99;
los_pv2 = zeros(1,size(spv2,1));
bs = (max(spv2(:)) - min(spv2(:))) * mlt + min(spv2(:));
bs = (prctile(spv2(:),p2) - prctile(spv2(:),p1)) * mlt + prctile(spv2(:),p1);
for k = 1:size(spv2,1)
    los_pv2(k) = sum(spv2(k,:)<bs);
end

los_som2 = zeros(1,size(ssom2,1));
bs = (max(ssom2(:)) - min(ssom2(:))) * mlt + min(ssom2(:));
bs = (prctile(ssom2(:),p2) - prctile(ssom2(:),p1)) * mlt + prctile(ssom2(:),p1);
for k = 1:size(ssom2,1)
    los_som2(k) = sum(ssom2(k,:)<bs);
end

los_pvn2 = zeros(1,size(pvn2,1));
bs = (max(pvn2(:)) - min(pvn2(:))) * mlt + min(pvn2(:));
bs = (prctile(pvn2(:),p2) - prctile(pvn2(:),p1)) * mlt + prctile(pvn2(:),p1);
for k = 1:size(pvn2,1)
    los_pvn2(k) = sum(pvn2(k,:)<bs);
end

los_somn2 = zeros(1,size(somn2,1));
bs = (max(somn2(:)) - min(somn2(:))) * mlt + min(somn2(:));
bs = (prctile(somn2(:),p2) - prctile(somn2(:),p1)) * mlt + prctile(somn2(:),p1);
for k = 1:size(somn2,1)
    los_somn2(k) = sum(somn2(k,:)<bs);
end

los_nt2 = zeros(1,size(snt2,1));
bs = (max(snt2(:)) - min(snt2(:))) * mlt + min(snt2(:));
bs = (prctile(snt2(:),p2) - prctile(snt2(:),p1)) * mlt + prctile(snt2(:),p1);
for k = 1:size(snt2,1)
    los_nt2(k) = sum(snt2(k,:)<bs);
end

figure
plot(ntd,los_nt2,'.','MarkerSize',20,'Color',[0.7 0.7 0.7])
hold on
plot(somnd,los_somn2,'c.','MarkerSize',20)
plot(spvd,los_pv2,'r.','MarkerSize',20)
plot(ssomd,los_som2,'b.','MarkerSize',20)

%% MDA

v1 = [spvd1 ssomd1 sntd1]';
v1 = (v1 - mean(v1)) / std(v1);
v1b = [spvd2 ssomd2 sntd2]';
v1b = (v1b - mean(v1b)) / std(v1b);
v2 = [spve1 ssome1 snte1]';
v2 = (v2 - mean(v2)) / std(v2);
v2b = [spve2 ssome2 snte2]';
v2b = (v2b - mean(v2b)) / std(v2b);
v3 = [si_pv2 si_som2 si_snt2]';
v3 = (v3 - mean(v3)) / std(v3);
v4 = [los_pv1 los_som1 los_nt1]';
v4 = (v4 - mean(v4)) / std(v4);
v5 = [los_pv2 los_som2 los_nt2]';
v5 = (v5 - mean(v5)) / std(v5);
vars = [v1 v2 v1b v2b...   % corr. w PV, corr. w SOM.
    v3...    % supression index (HZIn)
    v4 v5];  % supression time (HZOut), supression time (HZIn)

celltypes = {1:sizepv(1); sizepv(1)+1:sizepv(1)+sizesom(1); sizepv(1)+sizesom(1)+1:size(vars,1)};
celltyps = [zeros(sizepv(1),1); ones(sizesom(1),1); ones(size(snt1,1),1)*2];

W = lda2(vars,celltyps);
scores = [ones(length(celltyps),1) vars] * W';
score_1 = scores(:,2);
score_2 = scores(:,3);

% if PERM == 1
%     score = Xd(indREV,:);
% else
%     score = Xd;
% end

figure
plot(score_1(celltypes{3}),score_2(celltypes{3}),'.','MarkerSize',20,'Color',[0.7 0.7 0.7]);
hold on
plot(score_1(celltypes{1}),score_2(celltypes{1}),'r.','MarkerSize',20);
plot(score_1(celltypes{2}),score_2(celltypes{2}),'b.','MarkerSize',20);

[coeff PCAscores] = princomp([score_1 score_2]);
PC1 = PCAscores(:,1);
PC2 = PCAscores(:,2);

figure
plot(PC1(celltypes{3}),PC2(celltypes{3}),'.','MarkerSize',20,'Color',[0.7 0.7 0.7]);
hold on
plot(PC1(celltypes{1}),PC2(celltypes{1}),'r.','MarkerSize',20);
plot(PC1(celltypes{2}),PC2(celltypes{2}),'b.','MarkerSize',20);

%% MDA - PV pairs and SOM pairs included as two groups

v1 = [spvd1 ssomd1 pvnd1 somnd1]';
v1 = (v1 - mean(v1)) / std(v1);
v2 = [spve1 ssome1 pvne1 somne1]';
v2 = (v2 - mean(v2)) / std(v2);
v3 = [si_pv2 si_som2 si_pvn2 si_somn2]';
v3 = (v3 - mean(v3)) / std(v3);
v4 = [los_pv1 los_som1 los_pvn1 los_somn1]';
v4 = (v4 - mean(v4)) / std(v4);
v5 = [los_pv2 los_som2 los_pvn2 los_somn2]';
v5 = (v5 - mean(v5)) / std(v5);
vars = [v1 v2...   % corr. w PV, corr. w SOM.
    v3...    % supression index (HZIn)
    v4 v5];  % supression time (HZOut), supression time (HZIn)
celltyps = [zeros(sizepv(1),1); ones(sizesom(1),1); ones(size(pvn1,1),1)*2; ones(size(somn1,1),1)*3];

W = lda2(vars,celltyps);
v1 = [spvd1 ssomd1 sntd1]';
v1 = (v1 - mean(v1)) / std(v1);
v2 = [spve1 ssome1 snte1]';
v2 = (v2 - mean(v2)) / std(v2);
v3 = [si_pv2 si_som2 si_snt2]';
v3 = (v3 - mean(v3)) / std(v3);
v4 = [los_pv1 los_som1 los_nt1]';
v4 = (v4 - mean(v4)) / std(v4);
v5 = [los_pv2 los_som2 los_nt2]';
v5 = (v5 - mean(v5)) / std(v5);
vars2 = [v1 v2...   % corr. w PV, corr. w SOM.
    v3...    % supression index (HZIn)
    v4 v5];  % supression time (HZOut), supression time (HZIn)

scores = [ones(size(vars2,1),1) vars2] * W';
score_1 = scores(:,2);
score_2 = scores(:,3);

celltypes = {1:sizepv(1); sizepv(1)+1:sizepv(1)+sizesom(1); sizepv(1)+sizesom(1)+1:size(vars2,1)};

% if PERM == 1
%     score = Xd(indREV,:);
% else
%     score = Xd;
% end

figure
plot(score_1(celltypes{1}),score_2(celltypes{1}),'r.','MarkerSize',20);
hold on
plot(score_1(celltypes{2}),score_2(celltypes{2}),'b.','MarkerSize',20);
plot(score_1(celltypes{3}),score_2(celltypes{3}),'.','MarkerSize',20,'Color',[0.7 0.7 0.7]);

[coeff PCAscores] = princomp([score_1 score_2]);
PC1 = PCAscores(:,1);
PC2 = PCAscores(:,2);

figure
plot(PC1(celltypes{3}),PC2(celltypes{3}),'.','MarkerSize',20,'Color',[0.7 0.7 0.7]);
hold on
plot(PC1(celltypes{1}),PC2(celltypes{1}),'r.','MarkerSize',20);
plot(PC1(celltypes{2}),PC2(celltypes{2}),'b.','MarkerSize',20);

%% MDA - PV pairs and SOM pairs included as one group

v1 = [spvd1 ssomd1 pvnd1 somnd1]';
v1 = (v1 - mean(v1)) / std(v1);
v1b = [spvd2 ssomd2 pvnd2 somnd2]';
v1b = (v1b - mean(v1b)) / std(v1b);
v2 = [spve1 ssome1 pvne1 somne1]';
v2 = (v2 - mean(v2)) / std(v2);
v2b = [spve2 ssome2 pvne2 somne2]';
v2b = (v2b - mean(v2b)) / std(v2b);
v3 = [si_pv2 si_som2 si_pvn2 si_somn2]';
v3 = (v3 - mean(v3)) / std(v3);
v4 = [los_pv1 los_som1 los_pvn1 los_somn1]';
v4 = (v4 - mean(v4)) / std(v4);
v5 = [los_pv2 los_som2 los_pvn2 los_somn2]';
v5 = (v5 - mean(v5)) / std(v5);
vars = [v1 v2 v1b v2b...   % corr. w PV, corr. w SOM.
    v3...    % supression index (HZIn)
    v4 v5];  % supression time (HZOut), supression time (HZIn)
celltyps = [zeros(sizepv(1),1); ones(sizesom(1),1); ones(size(pvn1,1),1)*2; ones(size(somn1,1),1)*2];

W = lda2(vars,celltyps);
v1 = [spvd1 ssomd1 pvnd1 somnd1 sntd1]';
v1 = (v1 - mean(v1)) / std(v1);
v1b = [spvd2 ssomd2 pvnd2 somnd2 sntd2]';
v1b = (v1b - mean(v1b)) / std(v1b);
v2 = [spve1 ssome1 pvne1 somne1 snte1]';
v2 = (v2 - mean(v2)) / std(v2);
v2b = [spve2 ssome2 pvne2 somne2 snte2]';
v2b = (v2b - mean(v2b)) / std(v2b);
v3 = [si_pv2 si_som2 si_pvn2 si_somn2 si_snt2]';
v3 = (v3 - mean(v3)) / std(v3);
v4 = [los_pv1 los_som1 los_pvn1 los_somn1 los_nt1]';
v4 = (v4 - mean(v4)) / std(v4);
v5 = [los_pv2 los_som2 los_pvn2 los_somn2 los_nt2]';
v5 = (v5 - mean(v5)) / std(v5);
vars2 = [v1 v2 v1b v2b...   % corr. w PV, corr. w SOM.
    v3...    % supression index (HZIn)
    v4 v5];  % supression time (HZOut), supression time (HZIn)

scores = [ones(size(vars2,1),1) vars2] * W';
score_1 = scores(:,2);
score_2 = scores(:,3);

celltypes = {1:sizepv(1); sizepv(1)+1:sizepv(1)+sizesom(1); ...
    sizepv(1)+sizesom(1)+1:sizepv(1)+sizesom(1)+size(pvn1,1)+size(somn1,1); ...
    sizepv(1)+sizesom(1)+size(pvn1,1)+size(somn1,1)+1:size(vars2,1)};

% if PERM == 1
%     score = Xd(indREV,:);
% else
%     score = Xd;
% end

figure
plot(score_1(celltypes{4}),score_2(celltypes{4}),'.','MarkerSize',15,'Color',[0.7 0.7 0.7]);
hold on
plot(score_1(celltypes{3}),score_2(celltypes{3}),'k.','MarkerSize',20);
plot(score_1(celltypes{1}),score_2(celltypes{1}),'r.','MarkerSize',20);
plot(score_1(celltypes{2}),score_2(celltypes{2}),'b.','MarkerSize',20);

[coeff PCAscores] = princomp([score_1 score_2]);
PC1 = PCAscores(:,1);
PC2 = PCAscores(:,2);

figure
plot(PC1(celltypes{4}),PC2(celltypes{4}),'.','MarkerSize',15,'Color',[0.7 0.7 0.7]);
hold on
plot(PC1(celltypes{3}),PC2(celltypes{3}),'k.','MarkerSize',20);
plot(PC1(celltypes{1}),PC2(celltypes{1}),'r.','MarkerSize',20);
plot(PC1(celltypes{2}),PC2(celltypes{2}),'b.','MarkerSize',20);

%% temp

figure
plot(los_nt1,los_nt2,'.','MarkerSize',20,'Color',[0.7 0.7 0.7])
hold on
plot(los_pv1,los_pv2,'r.','MarkerSize',20)
plot(los_som1,los_som2,'b.','MarkerSize',20)