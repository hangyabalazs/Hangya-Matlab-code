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

% pv_MatrixPsth = [pv_MatrixPsth1 pv_MatrixPsth2];
% som_MatrixPsth = [som_MatrixPsth1 som_MatrixPsth2];
% non_tagged_MatrixPsth = [non_tagged_MatrixPsth1 non_tagged_MatrixPsth2];

%% remove NaNs

[x1 y1] = find(isnan(non_tagged_MatrixPsth1));
[x2 y2] = find(isnan(non_tagged_MatrixPsth2));
x = union(unique(x1),unique(x2));
non_tagged_MatrixPsth1(x,:) = [];
non_tagged_MatrixPsth2(x,:) = [];

%% restrict in time - new data

% pv_MatrixPsth1 = pv_MatrixPsth1(:,57:157);
% som_MatrixPsth1 = som_MatrixPsth1(:,57:157);
% non_tagged_MatrixPsth1 = non_tagged_MatrixPsth1(:,57:157);
% pv_MatrixPsth2 = pv_MatrixPsth2(:,57:157);
% som_MatrixPsth2 = som_MatrixPsth2(:,57:157);
% non_tagged_MatrixPsth2 = non_tagged_MatrixPsth2(:,57:157);
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

%% restrict in time

spv1 = spv1(:,57:157);
ssom1 = ssom1(:,57:157);
snt1 = snt1(:,57:157);
spv2 = spv2(:,57:157);
ssom2 = ssom2(:,57:157);
snt2 = snt2(:,57:157);

%% concatenate

spv = [spv1 spv2];
ssom = [ssom1 ssom2];
snt = [snt1 snt2];
pvn = [pvn1 pvn2];
somn = [somn1 somn2];

figure;plot(spv','r')
hold on;plot(ssom','b')
plot(median(spv),'r','LineWidth',3)
plot(median(ssom),'b','LineWidth',3)

%% PCA - all included

[coeff scores] = princomp([spv; ssom; pvn; somn]);
inx = [ones(sizepv(1),1); ones(sizesom(1),1)*2; ones(size(pvn,1),1)*3; ones(size(somn,1),1)*4];
for k = 1:4
    PC{k} = scores(:,k);
    pvPC{k} = PC{k}(inx==1);
    somPC{k} = PC{k}(inx==2);
    pvnPC{k} = PC{k}(inx==3);
    somnPC{k} = PC{k}(inx==4);
end

figure
plot(pvnPC{2},pvnPC{3},'.','MarkerSize',20,'Color','y')
hold on
plot(somnPC{2},somnPC{3},'.','MarkerSize',20,'Color','c')
plot(pvPC{2},pvPC{3},'r.','MarkerSize',20)
plot(somPC{2},somPC{3},'b.','MarkerSize',20)


%% distance from mean

medpv = mean(spv1);
spvd = zeros(1,sizepv(1));
for k = 1:sizepv(1)
    spvd(k) = sum((spv1(k,:)-medpv).^2);
end

pvnd = zeros(1,size(pvn,1));
for k = 1:size(pvn,1)
    pvnd(k) = sum((pvn1(k,:)-medpv).^2);
end

medpvn = mean(pvn1);
spve = zeros(1,sizepv(1));
for k = 1:sizepv(1)
    spve(k) = sum((spv1(k,:)-medpvn).^2);
end

pvne = zeros(1,size(pvn,1));
for k = 1:size(pvn,1)
    pvne(k) = sum((pvn1(k,:)-medpvn).^2);
end

figure
plot(spvd,spve,'r.','MarkerSize',20)
hold on
plot(pvnd,pvne,'y.','MarkerSize',20)

%% distance from mean

medpv = mean(spv1);
spvd = zeros(1,sizepv(1));
for k = 1:sizepv(1)
    spvd(k) = sum((spv1(k,:)-medpv).^2);
end

ssomd = zeros(1,sizesom(1));
for k = 1:sizesom(1)
    ssomd(k) = sum((ssom1(k,:)-medpv).^2);
end

pvnd = zeros(1,size(pvn,1));
for k = 1:size(pvn,1)
    pvnd(k) = sum((pvn1(k,:)-medpv).^2);
end

medsom = mean(ssom1);
spve = zeros(1,sizepv(1));
for k = 1:sizepv(1)
    spve(k) = sum((spv1(k,:)-medsom).^2);
end

ssome = zeros(1,sizesom(1));
for k = 1:sizesom(1)
    ssome(k) = sum((ssom1(k,:)-medsom).^2);
end

pvne = zeros(1,size(pvn,1));
for k = 1:size(pvn,1)
    pvne(k) = sum((pvn1(k,:)-medsom).^2);
end

figure
plot(spvd,spve,'r.','MarkerSize',20)
hold on
plot(ssomd,ssome,'b.','MarkerSize',20)
plot(pvnd,pvne,'y.','MarkerSize',20)

%% multiply with template

template = mean(spv1([1:7 11 12],:));
spvd = zeros(1,sizepv(1));
for k = 1:sizepv(1)
    spvd(k) = sum(spv1(k,:).*template);
end

pvnd = zeros(1,size(pvn,1));
for k = 1:size(pvn,1)
    pvnd(k) = sum(pvn1(k,:).*template);
end

ssomd = zeros(1,sizesom(1));
for k = 1:sizesom(1)
    ssomd(k) = sum(ssom1(k,:).*template);
end

somnd = zeros(1,size(somn,1));
for k = 1:size(somn,1)
    somnd(k) = sum(somn1(k,:).*template);
end

template2 = mean(pvn1(1:23,:));
spve = zeros(1,sizepv(1));
for k = 1:sizepv(1)
    spve(k) = sum(spv1(k,:).*template2);
end

pvne = zeros(1,size(pvn,1));
for k = 1:size(pvn,1)
    pvne(k) = sum(pvn1(k,:).*template2);
end

ssome = zeros(1,sizesom(1));
for k = 1:sizesom(1)
    ssome(k) = sum(ssom1(k,:).*template2);
end

somne = zeros(1,size(somn,1));
for k = 1:size(somn,1)
    somne(k) = sum(somn1(k,:).*template2);
end

figure
plot(spvd,spve,'r.','MarkerSize',20)
hold on
plot(pvnd,pvne,'y.','MarkerSize',20)
plot(ssomd,ssome,'b.','MarkerSize',20)
plot(somnd,somne,'c.','MarkerSize',20)

%% Voronoi diagram

x = [[spvd'; ssomd'; -100; -100; 100; 100] [spve'; ssome'; -100; 100; -100; 100]];
[v,c] = voronoin(x); 
for pg = 1:length(c)
    if all(c{pg}~=1)   % If at least one of the indices is 1, then it is an open region and we can't patch that.
        ip = inpolygon(x(:,1),x(:,2),v(c{pg},1),v(c{pg},2));
        if ismember(find(ip),1:sizepv(1))
            edgeclr = [1 0 0];
            faceclr = [0.75 0.25 0.25];
        elseif ismember(find(ip),sizepv(1)+1:sizepv(1)+sizesom(1))
            edgeclr = [0 0 1];
            faceclr = [0.25 0.25 0.75];
        end
        patch(v(c{pg},1),v(c{pg},2),edgeclr,'FaceColor',faceclr,'FaceAlpha',0.5); % use color i.
    end
end

%% multiply with template

template = mean(spv1([1:7 11 12],:));
spvd2 = zeros(1,sizepv(1));
for k = 1:sizepv(1)
    spvd2(k) = sum(spv2(k,:).*template);
end

pvnd2 = zeros(1,size(pvn,1));
for k = 1:size(pvn,1)
    pvnd2(k) = sum(pvn2(k,:).*template);
end

ssomd2 = zeros(1,sizesom(1));
for k = 1:sizesom(1)
    ssomd2(k) = sum(ssom2(k,:).*template);
end

somnd2 = zeros(1,size(somn,1));
for k = 1:size(somn,1)
    somnd2(k) = sum(somn2(k,:).*template);
end

template2 = mean(pvn1(1:23,:));
spve2 = zeros(1,sizepv(1));
for k = 1:sizepv(1)
    spve2(k) = sum(spv2(k,:).*template2);
end

pvne2 = zeros(1,size(pvn,1));
for k = 1:size(pvn,1)
    pvne2(k) = sum(pvn2(k,:).*template2);
end

ssome2 = zeros(1,sizesom(1));
for k = 1:sizesom(1)
    ssome2(k) = sum(ssom2(k,:).*template2);
end

somne2 = zeros(1,size(somn,1));
for k = 1:size(somn,1)
    somne2(k) = sum(somn2(k,:).*template2);
end

figure
plot(spvd2,spve2,'r.','MarkerSize',20)
hold on
plot(pvnd2,pvne2,'y.','MarkerSize',20)
plot(ssomd2,ssome2,'b.','MarkerSize',20)
plot(somnd2,somne2,'c.','MarkerSize',20)

%% template approach, include both events

figure
plot(spvd+spve2,spve+spvd2,'r.','MarkerSize',20)
hold on
plot(pvnd+pvne2,pvne+pvnd2,'y.','MarkerSize',20)
plot(ssomd+ssome2,ssome+ssomd2,'b.','MarkerSize',20)
plot(somnd+somne2,somne+somnd2,'c.','MarkerSize',20)

%% firing in mid and side portions

mid_pv = sum(spv1(:,51:125),2);
mid_som = sum(ssom1(:,51:125),2);
mid_pvn = sum(pvn1(:,51:125),2);

side_pv = sum(spv1(:,[1:50 126:201]),2);
side_som = sum(ssom1(:,[1:50 126:201]),2);
side_pvn = sum(pvn1(:,[1:50 126:201]),2);

figure
plot(mid_pv,side_pv,'r.','MarkerSize',20)
hold on
plot(mid_pvn,side_pvn,'y.','MarkerSize',20)
plot(mid_som,side_som,'b.','MarkerSize',20)

%% bigger > peak-valley/2 (pv-pair: peak; som: valley)

mid_pv = sum(spv1(:,51:125),2);
mid_som = sum(ssom1(:,51:125),2);
mid_pvn = sum(pvn1(:,51:125),2);

mlt = 0.5;
pov_pv = zeros(1,size(spv1,1));
for k = 1:size(spv1,1)
    bs = (max(spv1(k,:)) - min(spv1(k,:))) * mlt + min(spv1(k,:));
    pov_pv(k) = sum(spv1(k,:)>bs);
end

pov_som = zeros(1,size(ssom1,1));
for k = 1:size(ssom1,1)
    bs = (max(ssom1(k,:)) - min(ssom1(k,:))) * mlt + min(ssom1(k,:));
    pov_som(k) = sum(ssom1(k,:)>bs);
end

pov_pvn = zeros(1,size(pvn1,1));
for k = 1:size(pvn1,1)
    bs = (max(pvn1(k,:)) - min(pvn1(k,:))) * mlt + min(pvn1(k,:));
    pov_pvn(k) = sum(pvn1(k,:)>bs);
end

figure
plot(mid_pv,pov_pv,'r.','MarkerSize',20)
hold on
plot(mid_pvn,pov_pvn,'y.','MarkerSize',20)
plot(mid_som,pov_som,'b.','MarkerSize',20)

%% pov vs 

template = mean(spv1([1:2 4:7 11 12],:));
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

mlt = 0.2;
pov_pv = zeros(1,size(spv1,1));
for k = 1:size(spv1,1)
    bs = (max(spv1(k,:)) - min(spv1(k,:))) * mlt + min(spv1(k,:));
    pov_pv(k) = sum(spv1(k,:)>bs);
end

pov_som = zeros(1,size(ssom1,1));
for k = 1:size(ssom1,1)
    bs = (max(ssom1(k,:)) - min(ssom1(k,:))) * mlt + min(ssom1(k,:));
    pov_som(k) = sum(ssom1(k,:)>bs);
end

pov_pvn = zeros(1,size(pvn1,1));
for k = 1:size(pvn1,1)
    bs = (max(pvn1(k,:)) - min(pvn1(k,:))) * mlt + min(pvn1(k,:));
    pov_pvn(k) = sum(pvn1(k,:)>bs);
end

figure
plot(spvd,pov_pv,'r.','MarkerSize',20)
hold on
plot(ssomd,pov_som,'b.','MarkerSize',20)
plot(pvnd,pov_pvn,'y.','MarkerSize',20)

%% pov vs mod 

template = mean(spv1([1:2 4:7 11 12],:));
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

mlt = 0.5;
pov_pv = zeros(1,size(spv1,1));
for k = 1:size(spv1,1)
    bs = (max(spv1(k,:)) - min(spv1(k,:))) * mlt + min(spv1(k,:));
    pov_pv(k) = mean(spv1(k,spv1(k,:)>bs));
end

pov_som = zeros(1,size(ssom1,1));
for k = 1:size(ssom1,1)
    bs = (max(ssom1(k,:)) - min(ssom1(k,:))) * mlt + min(ssom1(k,:));
    pov_som(k) = mean(ssom1(k,ssom1(k,:)>bs));
end

pov_pvn = zeros(1,size(pvn1,1));
for k = 1:size(pvn1,1)
    bs = (max(pvn1(k,:)) - min(pvn1(k,:))) * mlt + min(pvn1(k,:));
    pov_pvn(k) = mean(pvn1(k,pvn1(k,:)>bs));
end

figure
plot(spvd,pov_pv,'r.','MarkerSize',20)
hold on
plot(ssomd,pov_som,'b.','MarkerSize',20)
plot(pvnd,pov_pvn,'y.','MarkerSize',20)

%% value at max grad

template = mean(spv1([1:2 4:7 11 12],:));
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

vmg_pv = zeros(1,size(spv1,1));
for k = 1:size(spv1,1)
    inx = diff(spv1(k,:)) == max(diff(spv1(k,:)));
    vmg_pv(k) = spv1(k,inx);
end

vmg_pvn = zeros(1,size(pvn1,1));
for k = 1:size(pvn1,1)
    inx = diff(pvn1(k,:)) == max(diff(pvn1(k,:)));
    vmg_pvn(k) = pvn1(k,inx);
end

vmg_som = zeros(1,size(ssom1,1));
for k = 1:size(ssom1,1)
    inx = diff(ssom1(k,:)) == max(diff(ssom1(k,:)));
    vmg_som(k) = ssom1(k,inx);
end

figure
plot(spvd,vmg_pv,'r.','MarkerSize',20)
hold on
plot(ssomd,vmg_som,'b.','MarkerSize',20)
plot(pvnd,vmg_pvn,'y.','MarkerSize',20)

%% sum of silence

template = mean(spv1([1:2 4:7 11 12],:));
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

mlt = 0.25;
los_pv = zeros(1,size(spv1,1));
bs = (max(spv1(:)) - min(spv1(:))) * mlt + min(spv1(:));
for k = 1:size(spv1,1)
    los_pv(k) = sum(spv1(k,spv1(k,:)<bs));
end

los_som = zeros(1,size(ssom1,1));
bs = (max(ssom1(:)) - min(ssom1(:))) * mlt + min(ssom1(:));
for k = 1:size(ssom1,1)
    los_som(k) = sum(ssom1(k,ssom1(k,:)<bs));
end

los_pvn = zeros(1,size(pvn1,1));
bs = (max(pvn1(:)) - min(pvn1(:))) * mlt + min(pvn1(:));
for k = 1:size(pvn1,1)
    los_pvn(k) = sum(pvn1(k,pvn1(k,:)<bs));
end

figure
plot(spvd,los_pv,'r.','MarkerSize',20)
hold on
plot(ssomd,los_som,'b.','MarkerSize',20)
plot(pvnd,los_pvn,'y.','MarkerSize',20)

%% length of silence - THIS IS GOOD

% template = mean(spv1([1:2 4:7 11 12],:));
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
%% multiply with template

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

template2 = 2*mean(spv2) - mean(ssom2);
spve2 = zeros(1,sizepv(1));
for k = 1:sizepv(1)
    spve2(k) = sum(spv2(k,:).*template2);
end

pvne2 = zeros(1,size(pvn2,1));
for k = 1:size(pvn2,1)
    pvne2(k) = sum(pvn2(k,:).*template2);
end

ssome2 = zeros(1,sizesom(1));
for k = 1:sizesom(1)
    ssome2(k) = sum(ssom2(k,:).*template2);
end

figure
plot(spvd2,spve2,'r.','MarkerSize',20)
hold on
plot(pvnd2,pvne2,'y.','MarkerSize',20)
plot(ssomd2,ssome2,'b.','MarkerSize',20)

%% LDA

% dimensions: correlation with mean pv, correlation with mean som
% looking for best separating dimension

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

inx1 = 125:150;
inx2 = 1:100;
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

figure
plot(spvd2,si_pv2,'r.','MarkerSize',20)
hold on
plot(pvnd2,si_pvn2,'y.','MarkerSize',20)
plot(ssomd2,si_som2,'b.','MarkerSize',20)

figure
[nm xout] = hist(si_som2);
stairs(xout,nm,'Color','b','LineWidth',3)
hold on
[nm xout] = hist(si_pv2);
stairs(xout,nm,'Color','r')

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