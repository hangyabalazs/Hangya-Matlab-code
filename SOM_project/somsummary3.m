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

%% length of silence - THIS IS GOOD

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

mlt = 0.25;   % 0.25: maybe best; 0.3 also good
los_pv = zeros(1,size(spv1,1));
bs = (max(spv1(:)) - min(spv1(:))) * mlt + min(spv1(:));
for k = 1:size(spv1,1)
    los_pv(k) = sum(spv1(k,:)<bs);
end

los_som = zeros(1,size(ssom1,1));
bs = (max(ssom1(:)) - min(ssom1(:))) * mlt + min(ssom1(:));
for k = 1:size(ssom1,1)
    los_som(k) = sum(ssom1(k,:)<bs);
end

los_pvn = zeros(1,size(pvn1,1));
bs = (max(pvn1(:)) - min(pvn1(:))) * mlt + min(pvn1(:));
for k = 1:size(pvn1,1)
    los_pvn(k) = sum(pvn1(k,:)<bs);
end

figure
plot(spvd,los_pv,'r.','MarkerSize',20)
hold on
plot(ssomd,los_som,'b.','MarkerSize',20)
plot(pvnd,los_pvn,'y.','MarkerSize',20)

%% ------------------------------------------------------------------------
% SECOND EVENT
% -------------------------------------------------------------------------
%% length of silence - GOOD Y AXIS

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
los_pv = zeros(1,size(spv2,1));
bs = (max(spv2(:)) - min(spv2(:))) * mlt + min(spv2(:));
for k = 1:size(spv2,1)
    los_pv(k) = sum(spv2(k,:)<bs);
end

los_som = zeros(1,size(ssom2,1));
bs = (max(ssom2(:)) - min(ssom2(:))) * mlt + min(ssom2(:));
for k = 1:size(ssom2,1)
    los_som(k) = sum(ssom2(k,:)<bs);
end

los_somn = zeros(1,size(somn2,1));
bs = (max(somn2(:)) - min(somn2(:))) * mlt + min(somn2(:));
for k = 1:size(somn2,1)
    los_somn(k) = sum(somn2(k,:)<bs);
end

los_pvn = zeros(1,size(pvn2,1));
bs = (max(pvn2(:)) - min(pvn2(:))) * mlt + min(pvn2(:));
for k = 1:size(pvn2,1)
    los_pvn(k) = sum(pvn2(k,:)<bs);
end

figure
plot(spvd,los_pv,'r.','MarkerSize',20)
hold on
plot(ssomd,los_som,'b.','MarkerSize',20)
plot(somnd,los_somn,'c.','MarkerSize',20)
% plot(pvnd,los_pvn,'y.','MarkerSize',20)

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

template2 = mean(ssom1);
spve1 = zeros(1,sizepv(1));
for k = 1:sizepv(1)
    spve1(k) = sum(spv1(k,:).*template2);
end

ssome1 = zeros(1,sizesom(1));
for k = 1:sizesom(1)
    ssome1(k) = sum(ssom1(k,:).*template2);
end

figure
plot(spvd1,spve1,'r.','MarkerSize',20)
hold on
plot(ssomd1,ssome1,'b.','MarkerSize',20)

template = mean(spv2);
spvd2 = zeros(1,sizepv(1));
for k = 1:sizepv(1)
    spvd2(k) = sum(spv2(k,:).*template);
end

ssomd2 = zeros(1,sizesom(1));
for k = 1:sizesom(1)
    ssomd2(k) = sum(ssom2(k,:).*template);
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

figure
plot(spvd2,spve2,'r.','MarkerSize',20)
hold on
plot(ssomd2,ssome2,'b.','MarkerSize',20)

%% PCA

[coeff scores] = princomp([[spvd1 ssomd1]' [spve1 ssome1]']);
PC1_1 = scores(:,1);
PC2_1 = scores(:,2);
coeff1_1 = coeff(:,1);
coeff1_2 = coeff(:,2);

figure
plot(PC1_1(1:sizepv(1)),PC2_1(1:sizepv(1)),'r.','MarkerSize',20)
hold on
plot(PC1_1(sizepv(1)+1:end),PC2_1(sizepv(1)+1:end),'b.','MarkerSize',20)

%% PCA

[coeff scores] = princomp([[spvd2 ssomd2]' [spve2 ssome2]']);
PC1_2 = scores(:,1);
PC2_2 = scores(:,2);
coeff2_1 = coeff(:,1);
coeff2_2 = coeff(:,2);

figure
plot(PC2_2(1:sizepv(1)),PC1_2(1:sizepv(1)),'r.','MarkerSize',20)
hold on
plot(PC2_2(sizepv(1)+1:end),PC1_2(sizepv(1)+1:end),'b.','MarkerSize',20)

%% PCA vs supression time

template = coeff1_1(1) * mean(spv1) + coeff1_1(2) * mean(ssom1);
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

mlt = 0.25;   % 0.25: maybe best; 0.3 also good
los_pv = zeros(1,size(spv1,1));
bs = (max(spv1(:)) - min(spv1(:))) * mlt + min(spv1(:));
for k = 1:size(spv1,1)
    los_pv(k) = sum(spv1(k,:)<bs);
end

los_som = zeros(1,size(ssom1,1));
bs = (max(ssom1(:)) - min(ssom1(:))) * mlt + min(ssom1(:));
for k = 1:size(ssom1,1)
    los_som(k) = sum(ssom1(k,:)<bs);
end

los_pvn = zeros(1,size(pvn1,1));
bs = (max(pvn1(:)) - min(pvn1(:))) * mlt + min(pvn1(:));
for k = 1:size(pvn1,1)
    los_pvn(k) = sum(pvn1(k,:)<bs);
end

figure
plot(spvd,los_pv,'r.','MarkerSize',20)
hold on
plot(ssomd,los_som,'b.','MarkerSize',20)
plot(pvnd,los_pvn,'y.','MarkerSize',20)


%% PCA vs supression time2

template = coeff2_2(1) * mean(spv1) + coeff2_2(2) * mean(ssom1);
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
los_pv = zeros(1,size(spv2,1));
bs = (max(spv2(:)) - min(spv2(:))) * mlt + min(spv2(:));
for k = 1:size(spv2,1)
    los_pv(k) = sum(spv2(k,:)<bs);
end

los_som = zeros(1,size(ssom2,1));
bs = (max(ssom2(:)) - min(ssom2(:))) * mlt + min(ssom2(:));
for k = 1:size(ssom2,1)
    los_som(k) = sum(ssom2(k,:)<bs);
end

los_somn = zeros(1,size(somn2,1));
bs = (max(somn2(:)) - min(somn2(:))) * mlt + min(somn2(:));
for k = 1:size(somn2,1)
    los_somn(k) = sum(somn2(k,:)<bs);
end

los_pvn = zeros(1,size(pvn2,1));
bs = (max(pvn2(:)) - min(pvn2(:))) * mlt + min(pvn2(:));
for k = 1:size(pvn2,1)
    los_pvn(k) = sum(pvn2(k,:)<bs);
end

figure
plot(spvd,los_pv,'r.','MarkerSize',20)
hold on
plot(ssomd,los_som,'b.','MarkerSize',20)
plot(somnd,los_somn,'c.','MarkerSize',20)

%% MDA

vars = [[spvd2 ssomd2]' [spve2 ssome2]'];
celltypes = {1:sizepv(1); sizepv(1)+1:size(vars,1)};

[Xv, Xd, V, d] = lda(vars,celltypes);
score2 = Xd;

vars = [[spvd1 ssomd1]' [spve1 ssome1]'];
celltypes = {1:sizepv(1); sizepv(1)+1:size(vars,1)};

[Xv, Xd, V, d] = lda(vars,celltypes);
score1 = Xd;

% if PERM == 1
%     score = Xd(indREV,:);
% else
%     score = Xd;
% end

figure;
plot(score1(1:sizepv(1)),score2(1:sizepv(1)),'r.','MarkerSize',20);
hold on
plot(score1(sizepv(1)+1:end),score2(sizepv(1)+1:end),'b.','MarkerSize',20);

%% MDA vs supression time

template = coeff1_1(1) * mean(spv1) + coeff1_1(2) * mean(ssom1);
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

mlt = 0.25;   % 0.25: maybe best; 0.3 also good
los_pv = zeros(1,size(spv1,1));
bs = (max(spv1(:)) - min(spv1(:))) * mlt + min(spv1(:));
for k = 1:size(spv1,1)
    los_pv(k) = sum(spv1(k,:)<bs);
end

los_som = zeros(1,size(ssom1,1));
bs = (max(ssom1(:)) - min(ssom1(:))) * mlt + min(ssom1(:));
for k = 1:size(ssom1,1)
    los_som(k) = sum(ssom1(k,:)<bs);
end

los_pvn = zeros(1,size(pvn1,1));
bs = (max(pvn1(:)) - min(pvn1(:))) * mlt + min(pvn1(:));
for k = 1:size(pvn1,1)
    los_pvn(k) = sum(pvn1(k,:)<bs);
end

figure
plot(spvd,los_pv,'r.','MarkerSize',20)
hold on
plot(ssomd,los_som,'b.','MarkerSize',20)
plot(pvnd,los_pvn,'y.','MarkerSize',20)


%% PCA vs supression time2

template = coeff2_2(1) * mean(spv1) + coeff2_2(2) * mean(ssom1);
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
los_pv = zeros(1,size(spv2,1));
bs = (max(spv2(:)) - min(spv2(:))) * mlt + min(spv2(:));
for k = 1:size(spv2,1)
    los_pv(k) = sum(spv2(k,:)<bs);
end

los_som = zeros(1,size(ssom2,1));
bs = (max(ssom2(:)) - min(ssom2(:))) * mlt + min(ssom2(:));
for k = 1:size(ssom2,1)
    los_som(k) = sum(ssom2(k,:)<bs);
end

los_somn = zeros(1,size(somn2,1));
bs = (max(somn2(:)) - min(somn2(:))) * mlt + min(somn2(:));
for k = 1:size(somn2,1)
    los_somn(k) = sum(somn2(k,:)<bs);
end

los_pvn = zeros(1,size(pvn2,1));
bs = (max(pvn2(:)) - min(pvn2(:))) * mlt + min(pvn2(:));
for k = 1:size(pvn2,1)
    los_pvn(k) = sum(pvn2(k,:)<bs);
end

figure
plot(spvd,los_pv,'r.','MarkerSize',20)
hold on
plot(ssomd,los_som,'b.','MarkerSize',20)
plot(somnd,los_somn,'c.','MarkerSize',20)