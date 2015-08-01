function imisorun_circle
%IMISORUN_CIRCLE   Slow oscillation propagation.
%   IMISORUN_CIRCLE calculates the number of circles and reciprocal
%   connections and compares them to random.
%
%   See also IMISORUN.

% Directories
global DATADIR
global DATAPATH
pat = 'oitin1_wittner';
bi = [DATADIR 'human_SO\' pat '\' pat '.jpg'];
eg = '148a';
inpdir = [DATAPATH 'Ulbert\OITI_N1_EEG_' eg '\MImap\'];
try
    resdir = [DATAPATH 'Ulbert\OITI_N1_EEG_' eg '\Circle'];
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
nm_cols = 4;           % number of columns on the electrode grid
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

% Distance matrix
dm = zeros(chnum,chnum);
for x = 1:chnum
    for y = 1:chnum
        [i j] = gridind2sub(x,nm_cols);
        [k l] = gridind2sub(y,nm_cols);
        dm(x,y) = sqrt((i-k)^2 + (j-l)^2);
    end
end

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
nmax = 3;
for k = nmax:-1:1    % 20 most frequent
    r1 = rand(1,2);
    line(cG(trjc3_2(ico(k),:),1)'+[r1(1) rand(1,2) r1(1)]*10-5,...
        cG(trjc3_2(ico(k),:),2)'+[r1(2) rand(1,2) r1(2)]*10-5,...
        'Color',cm(round((nmax-k+1)/nmax*length(cm)),:),'LineWidth',2)
end
saveas(H,['mostfreqcirc3_' eg])
saveas(H,['mostfreqcirc3_' eg '.tif'])

mtn = 1000;
randcircdist = zeros(1,mtn);
for tno = 1:mtn
    [Vfrom_rand Vto_rand Vs_rand Vstart_rand Vend_rand ...
        Adj_rand rIMrand rIMLrand Nei_rand Neicell_rand] = ...
        randshuff(rIM,rIML,tmax,chnum,nm_arrows);
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
%     disp(['Number of 3-point circles: ' num2str(size(trjc3,1))])
%     disp(['Control: ' num2str(size(trjc3_rand,1))])
    randcircdist(tno) = size(trjc3_rand,1);
end
H = figure;
hist(randcircdist,20,'FaceColor',[ 0 0.4980 0],'EdgeColor',[ 0 0.4980 0])
hold on
line([size(trjc3,1) size(trjc3,1)],get(gca,'YLim'),'Color',[0.8471 0.1608 0],'LineWidth',4)
set(gcf,'Position',[440 372 560 138])
set(gca,'LineWidth',2)
set(gca,'YTick',[])
set(gca,'FontSize',20)
box off
saveas(H,['randcircdist_' eg])
saveas(H,['randcircdist_' eg '.tif'])

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
nmax = 10;
for k = nmax:-1:1    % 20 most frequent
    line(cG(trjc2_2(ico(k),1:2),1)'+rand(1,2)*10-5,...
        cG(trjc2_2(ico(k),1:2),2)'+rand(1,2)*10-5,...
        'Color',cm(round((nmax-k+1)/nmax*length(cm)),:),'LineWidth',2)
end
saveas(H,['mostfreqcirc2_' eg])
saveas(H,['mostfreqcirc2_' eg '.tif'])

mtn = 1000;
randcircdist = zeros(1,mtn);
for tno = 1:mtn
    [Vfrom_rand Vto_rand Vs_rand Vstart_rand Vend_rand ...
        Adj_rand rIMrand rIMLrand Nei_rand Neicell_rand] = ...
        randshuff(rIM,rIML,tmax,chnum,nm_arrows);
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
%     disp(['Number of reciprocal connections: ' num2str(size(trjc2,1))])
%     disp(['Control: ' num2str(size(trjc2_rand,1))])
    randrecdist(tno) = size(trjc2_rand,1);
end
H = figure;
hist(randrecdist,20,'FaceColor',[ 0 0.4980 0],'EdgeColor',[ 0 0.4980 0])
hold on
line([size(trjc2,1) size(trjc2,1)],get(gca,'YLim'),'Color',[0.8471 0.1608 0],'LineWidth',4)
set(gcf,'Position',[440 372 560 138])
set(gca,'LineWidth',2)
set(gca,'YTick',[])
set(gca,'FontSize',20)
box off
saveas(H,['randrecdist_' eg])
saveas(H,['randrecdist_' eg '.tif'])



% -------------------------------------------------------------------------
function [Vfrom_rand Vto_rand Vs_rand Vstart_rand Vend_rand ...
    Adj_rand rIMrand rIMLrand Nei_rand Neicell_rand] = ...
    randshuff(rIM,rIML,tmax,chnum,nm_arrows)

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

Nei_rand = zeros(chnum,chnum,nm_arrows);   % random electrode neighbourhood matrix
Neicell_rand = cell(1,nm_arrows);
for k = 1:nm_arrows
    ls1 = Vstart_rand > Vend_rand(k) & Vstart_rand < Vend_rand(k) + 100;    % next vector should start within 100 ms
    ls2 = Vfrom_rand == Vto_rand(k);
    inx = find(ls1&ls2);
    Nei_rand(Vfrom_rand(inx),Vto_rand(inx),k) = 1;
    Neicell_rand{k} = inx;
end



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

% -------------------------------------------------------------------------
function cG = gridloc_oiti40_mojzsis(chnum)

cG = zeros(chnum,2);
cG(1,:) = [621 264];
cG(2,:) = [652 323];
cG(3,:) = [676 375];
cG(4,:) = [703 432];
cG(5,:) = [726 481];
cG(6,:) = [558 296];
cG(7,:) = [584 341];
cG(8,:) = [612 393];
cG(9,:) = [639 452];
cG(10,:) = [668 506];
cG(11,:) = [506 323];
cG(12,:) = [533 371];
cG(13,:) = [563 418];
cG(14,:) = [587 477];
cG(15,:) = [621 526];
cG(16,:) = [446 334];
cG(17,:) = [474 389];
cG(18,:) = [505 441];
cG(19,:) = [535 500];
cG(20,:) = [564 555];