function imisorun_randconvdiv(pat,patno,eg,nm_rows,nm_cols)
%IMISORUN_RANDCONVDIV   Slow oscillation propagation.
%   IMISORUN_RANDCONVDIV calculates rondomized convergence and divergence
%   maps (see IMISORUN) restricted to 4x4 size for random convergence and
%   divergence strength distribution calculation.
%
%   See also IMISORUN, ICONVGROUPSTAT and ICALLER.

% Input argument check
error(nargchk(0,5,nargin))
if nargin < 5
    nm_cols = 5;           % number of columns on the electrode grid
end
if nargin < 4
    nm_rows = 4;           % number of rows on the electrode grid
end
if nargin < 3
    eg = '46';
end
if nargin < 2
    patno = num2str('40');
end
if nargin < 1
    pat = 'oiti40_mojzsis';
end

% Directories
global DATADIR
global DATAPATH
bi = [DATADIR 'human_SO\' pat '\' pat '.jpg'];
inpdir = [DATAPATH 'Ulbert\OITI_' patno '_EEG_' eg '\MImap\'];
try
    resdir = [DATAPATH 'Ulbert\OITI_' patno '_EEG_' eg '\Traject_randconvdiv'];
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
rIMr = rIM(1:16,1:16,:);
rIMLr = rIML(1:16,1:16,:);
rIMrand = zeros(size(rIMr));
rIMLrand = zeros(size(rIMLr));
for k = 1:tmax
    ak = rIMr(:,:,k);
    bk = rIMLr(:,:,k);
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
nm_arrows_rand = length(Vend_rand);

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
for k = 1:nm_arrows_rand
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
    if Vfrom(k) == mx && length(find((Vfrom==mx&((Vstart>=Vstart(k)&Vstart<=Vstart(k)+101)|(Vstart<=Vstart(k)&Vstart>=Vend(k)-101)))))>1
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