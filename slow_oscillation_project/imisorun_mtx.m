function imisorun_mtx
%IMISORUN_MTX   Slow oscillation propagation.
%   IMISORUN_MTX calculates
%   (i) preferred 3-step trajectories. Arrow start points should be within 
%       100 ms after the end point of the previous arrow.
%   (ii) origin of arrows. Arrow start points are considered if there is no
%       arrow income within 1 second
%   (iii) convergence of arrows. Arrows are considered if they arrive
%       within 100 ms.
%   (iv) divergence of arrows. Arrows are considered if they apart within
%       100 ms (i.e. in the same or in the next frame).
%   (v) arrow length, spreading time, velocity and coupling strenght
%       distribution.
%   Data matrices are saved.
%
%   See also IMISORUN.

% Directories
global DATADIR
global DATAPATH
pat = 'oiti39_gaal';
bi = [DATADIR 'human_SO\' pat '\' pat '.jpg'];
eg = '40';
inpdir = [DATAPATH 'Ulbert\OITI_39_EEG_' eg '\MImap\'];
try
    resdir = [DATAPATH 'Ulbert\OITI_39_EEG_' eg '\Mtx'];
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
nm_cols = 4;           % number of columns on the electrode grid (4 or 5)
if ~isequal(nm_rows*nm_cols,chnum)
    error('Improper number for grid electrodes.')
end
sr = 1000;             % sampling rate

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
FirstArrows = sorig1;
sorig2 = sum(orig2,3);
mxo2 = max(sorig2(:));
StartPoints = sorig2;

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
Trajectories = trj;

% 3-point Circles
trjc3 = trj(trj(:,1)==(trj(:,4)),:);
Circles3 = trjc3;

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
Circles2 = trjc2;

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
% sCnv = sum(Cnv,3);
Convergence = Cnv;

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
% sDiv = sum(Div,3);
Divergence = Div;

% Distributions
rIMLn = rIML(rIML>0);
rIMn = rIM(rIM>0);
ar = Adj .* repmat(dm,[1 1 tmax]);   % adjacency matrix with distances
arn = ar(ar>0);
dmm = dm / 100;     % in meters
rIMLs = rIML / sr;  % in seconds
V = repmat(dmm,[1 1 tmax]) ./ rIMLs;
V = V(~isnan(V)&~isinf(V)&(V>0));
Strength = rIMn;
Time = rIMLn;
Speed = V;
Distance = arn;

% Save
save(['mtxs_' eg],'FirstArrows','StartPoints','Trajectories','Circles2',...
    'Circles3','Convergence','Divergence','Strength','Time','Speed','Distance')

% -------------------------------------------------------------------------
function [inx1 inx2] = gridind2sub(ind,nm_cols)

inx1 = floor((ind-1)/nm_cols) + 1;
inx2 = rem(ind,nm_cols);
if inx2 == 0
    inx2 = nm_cols;
end