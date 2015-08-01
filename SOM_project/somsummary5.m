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

si_nt2 = zeros(1,size(snt2,1));
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

snte2 = zeros(1,size(snt2,1));
for k = 1:size(snt2,1)
    snte2(k) = sum(snt2(k,:).*template2);
end

figure
plot(spvd2,spve2,'r.','MarkerSize',20)
hold on
plot(ssomd2,ssome2,'b.','MarkerSize',20)
plot(sntd2,snte2,'k.','MarkerSize',20)

%% MDA

vars = [[spvd2 ssomd2 sntd2]' [spve2 ssome2 snte2]' [spvd1 ssomd1 sntd1]' [spve1 ssome1 snte1]'];
celltypes = {1:sizepv(1); sizepv(1)+1:sizepv(1)+sizesom(1); sizepv(1)+sizesom(1)+1:size(vars,1)};

[Xv, Xd, V, d] = lda(vars,celltypes);
score_1 = Xd;

% if PERM == 1
%     score = Xd(indREV,:);
% else
%     score = Xd;
% end

figure
[nm xout] = hist(score_1(celltypes{1}));
stairs(xout,nm/sum(nm),'Color','r','LineWidth',3)
hold on
[nm xout] = hist(score_1(celltypes{2}));
stairs(xout,nm/sum(nm),'Color','b')
[nm xout] = hist(score_1(celltypes{3}));
stairs(xout,nm/sum(nm),'Color',[0.7 0.7 0.7])

%% MDA

vars = [[spvd2 ssomd2 sntd2]' [spve2 ssome2 snte2]'];
celltypes = {1:sizepv(1); sizepv(1)+1:sizepv(1)+sizesom(1); sizepv(1)+sizesom(1)+1:size(vars,1)};

[Xv, Xd, V, d] = lda(vars,celltypes);
score2 = Xd;

vars = [[spvd1 ssomd1 sntd1]' [spve1 ssome1 snte1]'];
celltypes = {1:sizepv(1); sizepv(1)+1:sizepv(1)+sizesom(1); sizepv(1)+sizesom(1)+1:size(vars,1)};

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
plot(score1(sizepv(1)+1:sizepv(1)+sizesom(1)),score2(sizepv(1)+1:sizepv(1)+sizesom(1)),'b.','MarkerSize',20);
plot(score1(sizepv(1)+sizesom(1)+1:end),score2(sizepv(1)+sizesom(1)+1:end),'.','MarkerSize',20,'Color',[0.7 0.7 0.7]);

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

%% supression index

inx1 = 80:100;
inx2 = 81:120;

si_pv1 = zeros(1,sizepv(1));
for k = 1:sizepv(1)
    si_pv1(k) = sum(spv1(k,inx1));
end

si_som1 = zeros(1,sizesom(1));
for k = 1:sizesom(1)
    si_som1(k) = sum(ssom1(k,inx1));
end

si_nt1 = zeros(1,size(snt1,1));
for k = 1:size(snt1,1)
    si_nt1(k) = sum(snt1(k,inx1));
end

si_pv2 = zeros(1,sizepv(1));
for k = 1:sizepv(1)
    si_pv2(k) = sum(spv1(k,inx2));
end

si_som2 = zeros(1,sizesom(1));
for k = 1:sizesom(1)
    si_som2(k) = sum(ssom1(k,inx2));
end

si_nt2 = zeros(1,size(snt1,1));
for k = 1:size(snt1,1)
    si_nt2(k) = sum(snt1(k,inx2));
end

sj_pv1 = zeros(1,sizepv(1));
for k = 1:sizepv(1)
    sj_pv1(k) = sum(spv2(k,inx1));
end

sj_som1 = zeros(1,sizesom(1));
for k = 1:sizesom(1)
    sj_som1(k) = sum(ssom2(k,inx1));
end

sj_nt1 = zeros(1,size(snt2,1));
for k = 1:size(snt2,1)
    sj_nt1(k) = sum(snt2(k,inx1));
end

sj_pv2 = zeros(1,sizepv(1));
for k = 1:sizepv(1)
    sj_pv2(k) = sum(spv2(k,inx2));
end

sj_som2 = zeros(1,sizesom(1));
for k = 1:sizesom(1)
    sj_som2(k) = sum(ssom2(k,inx2));
end

sj_nt2 = zeros(1,size(snt2,1));
for k = 1:size(snt2,1)
    sj_nt2(k) = sum(snt2(k,inx2));
end

% figure
% plot(spvd2,si_pv2,'r.','MarkerSize',20)
% hold on
% plot(pvnd2,si_pvn2,'y.','MarkerSize',20)
% plot(ssomd2,si_som2,'b.','MarkerSize',20)
% 
% figure
% [nm xout] = hist(si_som2);
% stairs(xout,nm/sum(nm),'Color','b','LineWidth',3)
% hold on
% [nm xout] = hist(si_pv2);
% stairs(xout,nm/sum(nm),'Color','r')
% [nm xout] = hist(si_nt2);
% stairs(xout,nm/sum(nm),'Color',[0.7 0.7 0.7])

%% MDA

vars = [[si_pv2 si_som2 si_nt2]' [sj_pv2 sj_som2 sj_nt2]' ...
    [si_pv1 si_som1 si_nt1]' [sj_pv1 sj_som1 sj_nt1]'];
celltypes = {1:sizepv(1); sizepv(1)+1:sizepv(1)+sizesom(1); sizepv(1)+sizesom(1)+1:size(vars,1)};

[Xv, Xd, V, d] = lda(vars,celltypes);
score_2 = Xd;

% if PERM == 1
%     score = Xd(indREV,:);
% else
%     score = Xd;
% end

figure
[nm xout] = hist(score_2(celltypes{1}));
stairs(xout,nm/sum(nm),'Color','r','LineWidth',3)
hold on
[nm xout] = hist(score_2(celltypes{2}));
stairs(xout,nm/sum(nm),'Color','b')
[nm xout] = hist(score_2(celltypes{3}));
stairs(xout,nm/sum(nm),'Color',[0.7 0.7 0.7])

%% MDA vs MDA

figure
plot(score_1(celltypes{1}),score_2(celltypes{1}),'r.','MarkerSize',20);
hold on
plot(score_1(celltypes{2}),score_2(celltypes{2}),'b.','MarkerSize',20);
plot(score_1(celltypes{3}),score_2(celltypes{3}),'.','MarkerSize',20,'Color',[0.7 0.7 0.7]);