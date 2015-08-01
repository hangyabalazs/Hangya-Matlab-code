function imisorun_oiti37
%IMISORUN   Slow oscillation propagation.
%   IMISORUN calculates
%   (i) preferred 3-step trajectories. Arrow start points should be within 
%       100 ms after the end point of the previous arrow.
%   (ii) origin of arrows. Arrow start points are considered if there is no
%       arrow income within 1 second
%   (iii) convergence of arrows. Arrows are considered if they arrive
%       within 100 ms. Distribution of electrodes with converging input are
%       compared to random control with within-frame shuffling of arrows.
%   (iv) divergence of arrows. Arrows are considered if they apart within
%       100 ms (i.e. in the same or in the next frame). Distribution is
%       compared to random (see above).
%   (v) arrow length, spreading time, velocity and coupling strenght
%       distribution. Scatter plots showing the relation between the above
%       variables are also plotted.
%
%   See also IMISHIFT.

% Directories
global DATADIR
global DATAPATH
pat = 'oiti37_lukacs';
bi = [DATADIR 'human_SO\' pat '\' pat '.jpg'];
eg = '24';
inpdir = [DATAPATH 'Ulbert\OITI_37_EEG_' eg '\MImap\'];
try
    resdir = [DATAPATH 'Ulbert\OITI_37_EEG_' eg '\Traject'];
    cd(resdir)
catch
    mkdir(resdir)
    cd(resdir)
end

% Load
fn = [inpdir 'MIshiftfine.mat'];    % MI map
load(fn)

% Significance level
fn = [inpdir 'siglev_EEG' eg];
load(fn)
siglev = sl(4,2);   % sig. lev.: 0.0001

% Transform input variables
Adj = rIMax;        % adjacency matrix
Adj(rIMax<siglev) = NaN;
Adj(isnan(Adj)) = 0;
Adj(Adj>0) = 1;     % time-varying adjacency matrix
rIM = rIMax;
rIM(rIMax<siglev) = NaN;
rIM(isnan(rIM)) = 0;
rIML = rIMaxLoc;
rIML(rIMax<siglev) = NaN;
rIML(isnan(rIML)) = 0;

% Recording properties
chnum = size(Adj,1);   % number of grid electrodes
tmax = size(Adj,3);    % time axis length
nm_rows = 4;           % number of rows on the electrode grid
nm_cols = 5;           % number of columns on the electrode grid
if ~isequal(nm_rows*nm_cols,chnum)
    error('Improper number for grid electrodes.')
end
sr = 1000;             % sampling rate

% Loading the MR image
cG = eval(['gridloc_' pat '(chnum);']);    % electrode coordinates
cm = colormap;
I = imread(bi);         % read MR

% Arrow starts and ends
Astart = Adj .* permute(repmat(1:100:tmax*100,[chnum 1 chnum]),[1 3 2]);
Aend = Astart + rIML;
[Vfrom Vto Vs] = ind2sub(size(Adj),find(Adj));     % vector indeces
Vstart = (Vs - 1) * 100 + 1;
Vend = Vstart + rIML(find(Adj));
[xV yV zV] = ind2sub(size(Adj),find(Adj));     % vector indeces
nm_arrows = length(xV);    % number of vectors

% Random shuffling of arrows
rIMrand = zeros(size(rIM));
rIMLrand = zeros(size(rIML));
for k = 1:tmax
    ak = rIM(:,:,k);
    bk = rIML(:,:,k);
    rp = randperm(numel(ak));
    akr = zeros(size(ak));
    akr(:) = ak(rp(:));
    bkr = zeros(size(bk));
    bkr(:) = bk(rp(:));
    rIMrand(:,:,k) = akr;
    rIMLrand(:,:,k) = bkr;
end
[Vfrom_rand Vto_rand Vs_rand] = ind2sub(size(rIMrand),find(rIMrand));     % vector indeces
Vstart_rand = (Vs_rand - 1) * 100 + 1;
Vend_rand = Vstart_rand + rIMLrand(find(rIMrand));
Adj_rand = rIMrand;
Adj_rand(Adj_rand>0)=1;

% Electrode neighbourhood matrix
Nei = zeros(chnum,chnum,nm_arrows);
Neicell = cell(1,nm_arrows);
for k = 1:nm_arrows
    ls1 = Vstart > Vend(k) & Vstart < Vend(k) + 100;    % next vector should start within 100 ms
    ls2 = Vfrom == Vto(k);
    inx = find(ls1&ls2);
    Nei(Vfrom(inx),Vto(inx),k) = 1;
    Neicell{k} = inx;
end

Nei_rand = zeros(chnum,chnum,nm_arrows);   % random electrode neighbourhood matrix
Neicell_rand = cell(1,nm_arrows);
for k = 1:nm_arrows
    ls1 = Vstart_rand > Vend_rand(k) & Vstart_rand < Vend_rand(k) + 100;    % next vector should start within 100 ms
    ls2 = Vfrom_rand == Vto_rand(k);
    inx = find(ls1&ls2);
    Nei_rand(Vfrom_rand(inx),Vto_rand(inx),k) = 1;
    Neicell_rand{k} = inx;
end

% Distance matrix
dm = zeros(chnum,chnum);
for x = 1:chnum
    for y = 1:chnum
        [i j] = gridind2sub(x,nm_cols);
        [k l] = gridind2sub(y,nm_cols);
        dm(x,y) = sqrt((i-k)^2 + (j-l)^2);
    end
end
dm(18,:) = 0;
dm(:,18) = 0;

% Places of origin
orig1 = zeros(size(Adj));
orig2 = zeros(nm_rows,nm_cols,tmax);
for k = 1:nm_arrows
    if sum(sum(Adj(:,xV(k),max(zV(k)-4,1):zV(k)))) == 0
        orig1(xV(k),yV(k),zV(k)) = 1;
        [inx1 inx2] = gridind2sub(xV(k),nm_cols);
        orig2(inx1,inx2,zV(k)) = orig2(inx1,inx2,zV(k)) + 1;
    else
        orig1(xV(k),yV(k),zV(k)) = 0;
    end
end
sorig1 = sum(orig1,3);
mxo1 = max(sorig1(:));
H = figure;    % plot
imagesc(sorig1,[0 mxo1])       % 'first' arrows
saveas(H,['firstarrows_' eg])
sorig2 = sum(orig2,3);
mxo2 = max(sorig2(:));
H = figure;
imagesc(sorig2,[0 mxo2])       % starting points of 'first' arrows
saveas(H,['startpoints_' eg])

orig1_rand = zeros(size(Adj));   % random control
orig2_rand = zeros(nm_rows,nm_cols,tmax);
for k = 1:nm_arrows
    if sum(sum(Adj_rand(:,Vfrom_rand(k),max(Vs_rand(k)-4,1):Vs_rand(k)))) == 0
        orig1_rand(Vfrom_rand(k),Vto_rand(k),Vs_rand(k)) = 1;
        [inx1 inx2] = gridind2sub(Vfrom_rand(k),nm_cols);
        orig2_rand(inx1,inx2,Vs_rand(k)) = orig2(inx1,inx2,Vs_rand(k)) + 1;
    else
        orig1_rand(Vfrom_rand(k),Vto_rand(k),Vs_rand(k)) = 0;
    end
end
H = figure;    % plot
imagesc(sum(orig1_rand,3),[0 mxo1])       % 'first' arrows
saveas(H,['rnd_firstarrows_' eg])
H = figure;
imagesc(sum(orig2_rand,3),[0 mxo2])       % starting points of 'first' arrows
saveas(H,['rnd_startpoints_' eg])

% Trajectories
trj = [];
for k = 1:nm_arrows
    arr2 = Neicell{k};
    for kk = 1:length(arr2)
        arr3 = Neicell{arr2(kk)};
        lst = Vto(arr3);
        lenlst = length(lst);
        trjtemp = [ones(lenlst,1)*Vfrom(k) ones(lenlst,1)*Vto(k) ...
            ones(lenlst,1)*Vto(arr2(kk)) lst];
        trj = [trj; trjtemp];   % concatenating 3-step trajectories
    end
end

trj2 = [0 0 0 0];   % sorting different trajectories in descending order of freq.
countr = 0;
for k = 1:length(trj)
    trjdif = trj2 - repmat(trj(k,:),[size(trj2,1) 1]);
    inx = find(any(trjdif,2)==0);
    if ~isempty(inx)
        countr(inx) = countr(inx) + 1;
    else
        trj2 = [trj2; trj(k,:)];
        countr(end+1) = 1;
    end
end
[sco ico] = sort(countr,'descend');

H = figure;
imagesc(I)
hold on
next = 0;
k = 1;
nmax = 20;
while next < nmax
    if dm(trj2(ico(k),1),trj2(ico(k),4)) > 3    % 20 most frequent, >3 cm long trajectories
        next = next + 1;
        line(cG(trj2(ico(k),:),1)'+rand(1,4)*10-5,cG(trj2(ico(k),:),2)'+rand(1,4)*10-5,...
            'Color',cm(round((nmax-next+1)/nmax*length(cm)),:),'LineWidth',2)
    end
    k = k + 1;
end
saveas(H,['mostfreqtraj_' eg])

allxs = [];     % sum vector
allxe = [];
allys = [];
allye = [];
for k = 1:nm_arrows
    if dm(trj(k,1),trj(k,4)) > 3
        allxs = [allxs cG(trj(k,1),1)];
        allxe = [allxe cG(trj(k,4),1)];
        allys = [allys cG(trj(k,1),2)];
        allye = [allye cG(trj(k,4),2)];
    end
end
H = figure;
imagesc(I)
hold on
line([mean(allxs) mean(allxe)],[mean(allys) mean(allye)],'LineWidth',3)
plot(mean(allxs),mean(allys),'.','MarkerSize',20)
saveas(H,['sumvector_' eg])

H = figure;
imagesc(I)
hold on
for k = 1:min(size(trj,1),5000)     % plot all trajectories
    line([cG(trj(k,:),1)'+rand(1,4)*10],[cG(trj(k,:),2)'+rand(1,4)*10],...
        'Color',rand(1,3))
end
saveas(H,['alltraj_' eg])

% 3-point Circles
trjc3 = trj(trj(:,1)==(trj(:,4)),:);
trjc3_2 = [0 0 0 0];   % sorting different trajectories in descending order of freq.
countrc3 = 0;
for k = 1:length(trjc3)
    trjc3_dif = trjc3_2 - repmat(trjc3(k,:),[size(trjc3_2,1) 1]);
    inx = find(any(trjc3_dif,2)==0);
    if ~isempty(inx)
        countrc3(inx) = countrc3(inx) + 1;
    else
        trjc3_2 = [trjc3_2; trjc3(k,:)];
        countrc3(end+1) = 1;
    end
end
[sco ico] = sort(countrc3,'descend');

H = figure;
imagesc(I)
hold on
nmax = 20;
for k = nmax:-1:1    % 20 most frequent
    r1 = rand(1,2);
    line(cG(trjc3_2(ico(k),:),1)'+[r1(1) rand(1,2) r1(1)]*10-5,...
        cG(trjc3_2(ico(k),:),2)'+[r1(2) rand(1,2) r1(2)]*10-5,...
        'Color',cm(round((nmax-k+1)/nmax*length(cm)),:),'LineWidth',2)
end
saveas(H,['mostfreqcirc3_' eg])

H = figure;
imagesc(I)
hold on
for k = 1:size(trjc3,1)     % plot all circles
    line([cG(trjc3(k,:),1)'+rand(1,4)*10],[cG(trjc3(k,:),2)'+rand(1,4)*10],...
        'Color',rand(1,3))
end
saveas(H,['allcirc3_' eg])

trj_rand = [];       % random trajectories
for k = 1:nm_arrows
    arr2 = Neicell_rand{k};
    for kk = 1:length(arr2)
        arr3 = Neicell_rand{arr2(kk)};
        lst = Vto_rand(arr3);
        lenlst = length(lst);
        trjtemp = [ones(lenlst,1)*Vfrom_rand(k) ones(lenlst,1)*Vto_rand(k) ...
            ones(lenlst,1)*Vto_rand(arr2(kk)) lst];
        trj_rand = [trj_rand; trjtemp];   % concatenating 3-step trajectories
    end
end
trjc3_rand = trj_rand(trj_rand(:,1)==(trj_rand(:,4)),:);   % random 3-point circes
disp(['Number of 3-point circles: ' num2str(size(trjc3,1))])
disp(['Control: ' num2str(size(trjc3_rand,1))])

% Reciprocal connections
trjc2 = [];
for k = 1:nm_arrows
    arr2 = Neicell{k};
    for kk = 1:length(arr2)
        if Vfrom(k) == Vto(arr2(kk))
            trjc2_temp = [Vfrom(k) Vto(k) Vfrom(k)];
            trjc2 = [trjc2; trjc2_temp];
        end
    end
end

% trjc2 = [];
% for k = 1:nm_arrows
%     arr2 = Neicell{k};
%     lst = Vto(arr2);
%     lenlst = length(lst);
%     trjtemp = [ones(lenlst,1)*Vfrom(k) ones(lenlst,1)*Vto(k) lst];
%     trjc2 = [trjc2; trjtemp];
% end
% trjc2 = trjc2(trjc2(:,1)==(trjc2(:,3)),:);

trjc2_2 = [0 0 0];   % sorting different trajectories in descending order of freq.
countrc2 = 0;
for k = 1:length(trjc2)
    trjc2_dif = trjc2_2 - repmat(trjc2(k,:),[size(trjc2_2,1) 1]);
    inx = find(any(trjc2_dif,2)==0);
    if ~isempty(inx)
        countrc2(inx) = countrc2(inx) + 1;
    else
        trjc2_2 = [trjc2_2; trjc2(k,:)];
        countrc2(end+1) = 1;
    end
end
[sco ico] = sort(countrc2,'descend');

H = figure;
imagesc(I)
hold on
nmax = 20;
for k = nmax:-1:1    % 20 most frequent
    line(cG(trjc2_2(ico(k),1:2),1)'+rand(1,2)*10-5,...
        cG(trjc2_2(ico(k),1:2),2)'+rand(1,2)*10-5,...
        'Color',cm(round((nmax-k+1)/nmax*length(cm)),:),'LineWidth',2)
end
saveas(H,['mostfreqcirc2_' eg])

H = figure;
imagesc(I)
hold on
for k = 1:size(trjc2,1)     % plot all
    line([cG(trjc2(k,1:2),1)'+rand(1,2)*10-5],[cG(trjc2(k,1:2),2)'+rand(1,2)*10-5],...
        'Color',rand(1,3))
end
saveas(H,['allcirc2_' eg])

trjc2_rand = [];       % random
for k = 1:nm_arrows
    arr2 = Neicell_rand{k};
    for kk = 1:length(arr2)
        if Vfrom_rand(k) == Vto_rand(arr2(kk))
            trjc2_temp = [Vfrom_rand(k) Vto_rand(k) Vfrom_rand(k)];
            trjc2_rand = [trjc2_rand; trjc2_temp];
        end
    end
end
disp(['Number of reciprocal connections: ' num2str(size(trjc2,1))])
disp(['Control: ' num2str(size(trjc2_rand,1))])

% Convergence
tms = tmax * 100;   % arrows converging on an electrode within 100 ms
Cnv = zeros(nm_rows,nm_cols,tms);
for t = 1:tms
    for ltr = 1:chnum
        ls = Vto == ltr & Vend >= t & Vend < t + 100;
        fls = find(ls);
        [inx1 inx2] = gridind2sub(ltr,nm_cols);
        Cnv(inx1,inx2,t) = length(fls);
    end
end
H = figure;
sCnv = sum(Cnv,3);
mxsc = max(sCnv(:));
imagesc(sCnv,[0 mxsc])
saveas(H,['conv_' eg])

Cnv_rand = zeros(nm_rows,nm_cols,tms);   % random convergence
for t = 1:tms
    for ltr = 1:chnum
        ls = Vto_rand == ltr & Vend_rand >= t & Vend_rand < t + 100;
        fls = find(ls);
        [inx1 inx2] = gridind2sub(ltr,nm_cols);
        Cnv_rand(inx1,inx2,t) = length(fls);
    end
end
H = figure;
sCnv_rand = sum(Cnv_rand,3);
imagesc(sCnv_rand,[0 mxsc])
[nm xout] = hist(Cnv(:),(1:20));
[nm_rand xout_rand] = hist(Cnv_rand(:),(1:20));
saveas(H,['rnd_conv_' eg])
figure
bar(xout(2:end),nm(2:end))
hold on
plot(xout_rand(2:end),nm_rand(2:end),'k')

% arv = cell(1,chnum);      % distrution of income frequency on all electrodes
% for k = 1:nm_arrows
%     v = Vto(k);
%     arv{v}(end+1) = Vend(k);
% end
% add = [];
% for k = 1:chnum
%     add = [add diff(sort(arv{k}))];
% end
% figure
% hist(add(add<80),500)

mx = find(sCnv'==max(sCnv(:)));
H = figure;         % incoming arrows to the electrode with the largest no. of converging arrows
imagesc(I)
hold on
wf = zeros(1,chnum);
for k = 1:nm_arrows
    if Vto(k) == mx && any(Vto==mx&((Vend>Vend(k)&Vend<Vend(k)+100)|(Vend<Vend(k)&Vend>Vend(k)-100)))
        ltr = Vfrom(k);
        wf(ltr) = wf(ltr) + 1;
    end
end
[swf inx] = sort(wf);
for k0 = 1:chnum
    k = inx(k0);
    if wf(k) > 0
        line([cG(k,1) cG(mx,1)],[cG(k,2) cG(mx,2)],'LineWidth',3,'Color',...
            cm(round(wf(k)/max(wf)*length(cm)),:))
    end
end
saveas(H,['maxconv_' eg])

% Divergence
Div = zeros(nm_rows,nm_cols,tmax);      % arrows diverging from an electrodes in two condecutive frames
next = 1;
for t = 1:100:tmax*100
    for ltr = 1:chnum
        ls = Vfrom == ltr & Vstart >= t & Vstart <= t + 101;
        fls = find(ls);
        [inx1 inx2] = gridind2sub(ltr,nm_cols);
        Div(inx1,inx2,next) = length(fls);
    end
    next = next + 1;
end
H = figure;
sDiv = sum(Div,3);
imagesc(sDiv,[0 max(sDiv(:))])
saveas(H,['div_' eg])

Div_rand = zeros(nm_rows,nm_cols,tmax);     % random divergence
next = 1;
for t = 1:100:tmax*100
    for ltr = 1:chnum
        ls = Vfrom_rand == ltr & Vstart_rand >= t & Vstart_rand <= t + 101;
        fls = find(ls);
        [inx1 inx2] = gridind2sub(ltr,nm_cols);
        Div_rand(inx1,inx2,next) = length(fls);
    end
    next = next + 1;
end
H = figure;
sDiv_rand = sum(Div_rand,3);
imagesc(sDiv_rand,[0 max(sDiv(:))])
saveas(H,['rnd_div_' eg])
[nm xout] = hist(Div(:),(1:30));
[nm_rand xout_rand] = hist(Div_rand(:),(1:30));
figure
bar(xout(2:end),nm(2:end))
hold on
plot(xout_rand(2:end),nm_rand(2:end),'k')

mx = find(sDiv'==max(sDiv(:)));     % divergence from the elctrode with largest no. of diverging arrows
H = figure;
imagesc(I)
hold on
wt = zeros(1,chnum);
for k = 1:nm_arrows
    if Vfrom(k) == mx
        ltr = Vto(k);
        wt(ltr) = wt(ltr) + 1;
    end
end
cm = colormap;
[swt inx] = sort(wt);
for k0 = 1:chnum
    k = inx(k0);
    if wt(k) > 0
        line([cG(k,1) cG(mx,1)],[cG(k,2) cG(mx,2)],'LineWidth',3,'Color',...
            cm(max(round(wt(k)/max(wt)*length(cm)),1),:))
    end
end
saveas(H,['maxdiv_' eg])

% Distributions
H = figure;   % distribution of spreading time
rIMLn = rIML(rIML>0);
hist(rIMLn,100)
saveas(H,['timedist_' eg])

H = figure;   % distribution of MI(max)
rIMn = rIM(rIM>0);
hist(rIMn)
saveas(H,['MImaxdist_' eg])

H = figure;
ar = Adj .* repmat(dm,[1 1 tmax]);   % adjacency matrix with distances
arn = ar(ar>0);
[nm xout] = hist(arn,50);      % arrow length distribution
bar(xout,nm/sum(nm))
set(gca,'XTick',[1:5],'YTick',[0:0.1:1],'FontSize',20)
box off
ylim([0 0.4])
saveas(H,['lengthdist_' eg])
saveas(H,['lengthdist_' eg '.tif'])

H = figure;      % normalized arrow length distribution
[nm2 xout2] = hist(dm(dm>0),50);
if ~isequal(xout,xout2)
    error('Technical error 522.')
end
nmm = nm ./ nm2;    % normalization by the abundance of a given distance
nmm2 = nan2zero(nmm);
bar(xout,nmm2/sum(nmm2))
set(gca,'XTick',[1:5],'YTick',[0:0.1:1],'FontSize',20)
box off
saveas(H,['lengthdist_norm_' eg])
saveas(H,['lengthdist_norm_' eg '.tif'])

H = figure;     % frame no. vs. mean arrow length
arr = zeros(1,tmax);
for k = 1:tmax
    ark = ar(:,:,k);
    ark2 = ark(ark>0);
    arr(k) = mean(ark2);
end
plot(nan2zero(arr))

dmm = dm / 100;     % in meters
rIMLs = rIML / sr;  % in seconds
V = repmat(dmm,[1 1 tmax]) ./ rIMLs;
V = V(~isnan(V)&~isinf(V)&(V>0));
H = figure;
hist(V(:),1000)     % velocity
saveas(H,['velocitydist_' eg])

figure   % scatter plots
plot(rIMn,rIMLn,'.')
xlabel('Imax')
ylabel('max. loc.')

figure
plot(rIMLn,V,'.')
xlabel('max. loc.')
ylabel('velocity')

figure
plot(rIMn,V,'.')
xlabel('Imax')
ylabel('velocity')

figure
plot(arn,rIMLn,'.')
xlabel('arrow length')
ylabel('max. loc.')

figure
plot(arn,V,'.')
xlabel('arrow length')
ylabel('velocity')

figure
plot(arn,rIMn,'.')
xlabel('arrow length')
ylabel('Imax')



% -------------------------------------------------------------------------
function [inx1 inx2] = gridind2sub(ind,nm_cols)

inx1 = floor((ind-1)/nm_cols) + 1;
inx2 = rem(ind,nm_cols);
if inx2 == 0
    inx2 = nm_cols;
end

% -------------------------------------------------------------------------
function cG = gridloc_oiti37_lukacs(chnum)

cG = zeros(chnum,2);
cG(1,:) = [511 164];
cG(2,:) = [501 180];
cG(3,:) = [489 197];
cG(4,:) = [472 218];
cG(5,:) = [455 240];
cG(6,:) = [494 146];
cG(7,:) = [482 159];
cG(8,:) = [469 176];
cG(9,:) = [454 199];
cG(10,:) = [433 217];
cG(11,:) = [478 124];
cG(12,:) = [464 140];
cG(13,:) = [447 159];
cG(14,:) = [430 180];
cG(15,:) = [410 196];
cG(16,:) = [460 107];
cG(17,:) = [446 122];
cG(18,:) = [429 138];
cG(19,:) = [408 157];
cG(20,:) = [388 175];

% -------------------------------------------------------------------------
function cG = gridloc_oiti31_virag(chnum)

cG = zeros(chnum,2);
cG(1,:) = [517 166];
cG(2,:) = [509 157];
cG(3,:) = [497 147];
cG(4,:) = [480 130];
cG(5,:) = [496 182];
cG(6,:) = [488 173];
cG(7,:) = [476 160];
cG(8,:) = [461 143];
cG(9,:) = [480 195];
cG(10,:) = [473 188];
cG(11,:) = [460 176];
cG(12,:) = [444 160];
cG(13,:) = [463 211];
cG(14,:) = [456 203];
cG(15,:) = [444 192];
cG(16,:) = [429 176];

% -------------------------------------------------------------------------
function cG = gridloc_oiti39_gaal(chnum)

cG = zeros(chnum,2);
cG(1,:) = [385 154];
cG(2,:) = [401 141];
cG(3,:) = [416 126];
cG(4,:) = [427 114];
cG(5,:) = [364 140];
cG(6,:) = [380 126];
cG(7,:) = [395 112];
cG(8,:) = [410 100];
cG(9,:) = [341 127];
cG(10,:) = [357 112];
cG(11,:) = [377 99];
cG(12,:) = [387 87];
cG(13,:) = [318 115];
cG(14,:) = [332 101];
cG(15,:) = [348 86];
cG(16,:) = [364 73];

% -------------------------------------------------------------------------
function cG = gridloc_oitin1_wittner(chnum)

cG = zeros(chnum,2);
cG(1,:) = [529 498];
cG(2,:) = [498 516];
cG(3,:) = [481 539];
cG(4,:) = [461 564];
cG(5,:) = [501 467];
cG(6,:) = [467 481];
cG(7,:) = [445 506];
cG(8,:) = [424 528];
cG(9,:) = [469 428];
cG(10,:) = [437 446];
cG(11,:) = [410 469];
cG(12,:) = [387 496];
cG(13,:) = [442 399];
cG(14,:) = [414 420];
cG(15,:) = [388 442];
cG(16,:) = [359 466];