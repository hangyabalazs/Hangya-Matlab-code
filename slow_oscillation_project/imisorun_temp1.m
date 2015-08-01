function imisorun
%IMISORUN   Slow oscillation propagation.

% Load
fn = 'X:\In_Vivo\balazs\_analysis\Ulbert\EEG_12\MImap\MIshiftfine.mat';
load(fn)

% Transform input variables
rIM = rIMax;
rIML = rIMaxLoc;
Adj = rIMax;
siglev = 2.3475;   % sig. lev.: 0.0001
Adj(rIMax<siglev) = NaN;
rIML(rIMax<siglev) = NaN;
rIML(isnan(rIML)) = 0;
Adj(isnan(Adj)) = 0;
Tim = rIML;          
Adj(Adj>0) = 1;     % time-varying adjacency matrix

% Recording properties
chnum = size(Adj,1);   % number of grid electrodes
tmax = size(Adj,3);    % time axis length
nm_rows = 4;           % number of rows on the electrode grid
nm_cols = 5;           % number of columns on the electrode grid
if ~isequal(nm_rows*nm_cols,chnum)
    error('Improper number for grid electrodes.')
end
sr = 1000;             % sampling rate

% Arrow starts and ends
Astart = Adj .* permute(repmat(1:100:tmax*100,[chnum 1 chnum]),[1 3 2]);
Aend = Astart + Tim;

[Vfrom Vto Vs] = ind2sub(size(Adj),find(Adj));     % vector indeces
Vstart = (Vs - 1) * 100 + 1;
Vend = Vstart + rIML(find(Adj));


% Places of origin
[xV yV zV] = ind2sub(size(Adj),find(Adj));     % vector indeces
% orig1 = zeros(size(Adj));
% orig2 = zeros(nm_rows,nm_cols,tmax);
nm_arrows = length(xV);    % number of vectors
% for k = 1:nm_arrows
%     if sum(sum(Adj(:,xV(k),max(zV(k)-9,1):zV(k)))) == 0
%         orig1(xV(k),yV(k),zV(k)) = 1;
%         [inx1 inx2] = gridind2sub(xV(k),nm_cols);
%         orig2(inx1,inx2,zV(k)) = orig2(inx1,inx2,zV(k)) + 1;
%     else
%         orig1(xV(k),yV(k),zV(k)) = 0;
%     end
% end
% figure    % plot
% imagesc(sum(orig1,3))
% figure
% imagesc(sum(orig2,3))

% Trajectories
% for k = 1:tmax
%     frame = Adj(:,:,k);
%     nb1 = cell(1,chnum);
%     for ltr = 1:chnum
%         nb1{k} = find(frame(k,:))
%     end
% end

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

tic
trj = [];
vs = [];
for k = 1:nm_arrows
    arr2 = Neicell{k};
    for kk = 1:length(arr2)
        arr3 = Neicell{arr2(kk)};
        lst = Vto(arr3);
        lenlst = length(lst);
        trjtemp = [ones(lenlst,1)*Vfrom(k) ones(lenlst,1)*Vto(k) ...
            ones(lenlst,1)*Vto(arr2(kk)) lst];
        trj = [trj; trjtemp];
        vstemp = [k*ones(lenlst,1) arr2(kk)*ones(lenlst,1) arr3];
        vs = [vs; vstemp];
    end
end
toc
1

trj2 = [0 0 0 0];
countr = [0];
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
1





% -------------------------------------------------------------------------
function [inx1 inx2] = gridind2sub(ind,nm_cols)

inx1 = floor((ind-1)/nm_cols) + 1;
inx2 = rem(ind,nm_cols);
if inx2 == 0
    inx2 = 5;
end