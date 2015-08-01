function czplaceanalysis_Ntt
%CZPLACEANALYSIS_NTT   Analysis of spatial firing.
%   CZPLACEANALYSIS_NTT calculates place rate maps, autocorrelation,
%   crosscorrelation, place field similarity and complementarity index for
%   neuronal data. Edit code to modify input and output directories!
%
%   This version should be used for files derived from .Ntt using 
%   CZSPIKESORT!
%
%   See also CZPLACEANALYSIS and CZSPIKESORT.

% Directories & filename
global DATAPATH
fname = 'acin11s007_spk'
inpdir = [DATAPATH 'Czurko\discriminated2\new\' fname '\'];
resdir = [DATAPATH 'Czurko\discriminated2\new\' fname '\'];
fns = findstr(fname,'_');
titlestr = fname;
titlestr(fns) = ' ';
mm = pwd;
cd(resdir)
% load placecell_index
% load placemaps

% Import
fl = dir(inpdir);   % find text files
txts = {};
for k = 3:length(fl)
    if isequal(fl(k).name(end-2:end),'mat') && isequal(fl(k).name(1:2),'nr')
        txts{end+1} = fl(k).name;
    end
end
lentx = length(txts);

neuron = cell(1,lentx);     % load
for k = 1:lentx
    load([inpdir txts{k}])
    neuron{k} = data;
end
if lentx <= 9   % no. of rows and columns on place map plot
    rown = 3;
    coln = 3;
elseif lentx <= 20
    rown = 4;
    coln = 5;
elseif lentx <= 30
    rown = 5;
    coln = 6;
elseif lentx <= 42
    rown = 6;
    coln = 7;
else
    rown = 6;
    coln = 8;
    disp(fname)
    disp('Two many cells.')
end
nn = lentx;

% Read position data
load([DATAPATH 'Czurko\discriminated2\new\' fname '\x1_pos.mat']);
x_pos = data';
load([DATAPATH 'Czurko\discriminated2\new\' fname '\x1_pos_ts.mat']);
x_pos_ts = data';
load([DATAPATH 'Czurko\discriminated2\new\' fname '\y1_pos.mat']);
y_pos = data';
load([DATAPATH 'Czurko\discriminated2\new\' fname '\y1_pos_ts.mat']);
y_pos_ts = data';

% Place rate map
irhst = cell(1,nn);
rhst = cell(1,nn);
for k = 1:nn
    [irhst{k} rhst{k} thst] = czplace2(neuron{k},x_pos,y_pos,x_pos_ts);
    str=['subplot(rown,coln,' num2str(k) ');pcolor(nanpad(irhst{' num2str(k) '},1));shading flat'];     % plot
    eval(str)
    cl = get(gca,'CLim');
    set(gca,'CLim',[0 cl(2)])
    colorbar
    axis off
end
clear irhst
save placemaps rhst;      % save
set(gcf,'Position',[360 500 667 420])
title(titlestr)
saveas(gcf,'placemaps.fig')
saveas(gcf,'placemaps.jpg')
close(gcf)

% Place fields
pfsize = zeros(1,nn);
ms1 = cell(1,nn);
ms2 = cell(1,nn);
ms3 = cell(1,nn);
for k = 1:nn
    s = rhst{k};
    ms1{k} = s .* (zero2nan(double(s>b_mean_nonnan(s)+2*b_std_nonnan(s))));
    ms2{k} = s .* (zero2nan(double(s>b_mean_nonnan(s))));
    ms3{k} = s .* (zero2nan(double(s>b_max_nonnan(s)*0.5)));
    pfsize(k) = czpfsize(rhst{k});
end

% Selecting place cells (place field size >= 9 pixels) 
pci = find(pfsize>=9);
I = imread('placemaps.jpg');
figure
imagesc(I)
disp(['Place cells: ' num2str(pci)])
yn = [];
while ~any(strcmp(yn,{'y','n'}))
    yn = input('New place cells? [y/n]','s');
end
if isequal(yn,'y')
    pci = input('Place cells:');
end
placecell = neuron(pci);
np = length(pci);
save placecell_index pci    % save

% Autocorrelation
for k = 1:np
    czacorr(placecell{k})
    title(num2str(k))
    str = ['saveas(gcf,''autocorr' num2str(k) '.fig'')'];   % save
    eval(str)
    str = ['saveas(gcf,''autocorr' num2str(k) '.jpg'')'];
    eval(str)
end
close all

% Normalized cross-correlation
for x = 1:np
    for y = x+1:np
        [H1 H2 trsc] = czxcorr(placecell{x},placecell{y});   % window: +-50 ms
        title([num2str(x) '-' num2str(y)])
        str = ['saveas(gcf,''normcrosscorr_' num2str(x) '_' num2str(y) '.fig'')'];  % save
        eval(str)
        str = ['saveas(gcf,''normcrosscorr_' num2str(x) '_' num2str(y) '.jpg'')'];
        eval(str)
        figure(H1)
        title([num2str(x) '-' num2str(y)])
        str = ['saveas(gcf,''crosscorr_' num2str(x) '_' num2str(y) '.fig'')'];
        eval(str)
        str = ['saveas(gcf,''crosscorr_' num2str(x) '_' num2str(y) '.jpg'')'];
        eval(str)
        str = ['save(''trsc_' num2str(x) '_' num2str(y) ''',''trsc'')'];
        eval(str)
        close all
        
        [H1 H2] = czxcorr2(placecell{x},placecell{y});   % window: +-5 ms
        title([num2str(x) '-' num2str(y)])
        str = ['saveas(gcf,''nsmallxcorr_' num2str(x) '_' num2str(y) '.fig'')'];  % save
        eval(str)
        str = ['saveas(gcf,''nsmallxcorr_' num2str(x) '_' num2str(y) '.jpg'')'];
        eval(str)
        figure(H1)
        title([num2str(x) '-' num2str(y)])
        str = ['saveas(gcf,''smallxcorr_' num2str(x) '_' num2str(y) '.fig'')'];
        eval(str)
        str = ['saveas(gcf,''smallxcorr_' num2str(x) '_' num2str(y) '.jpg'')'];
        eval(str)
        close all
    end
end

% Linear correlation (place field similarity)
R = eye(np);
for x = 1:np
    for y = x+1:np
        R(x,y) = czpfs(rhst{pci(x)},rhst{pci(y)});
        R(y,x) = R(x,y);
    end
end
Rmod = eye(np);
for x = 1:np
    for y = x+1:np
        Rmod(x,y) = czpfs_mod(rhst{pci(x)},rhst{pci(y)});
        Rmod(y,x) = Rmod(x,y);
    end
end
save pfs R Rmod

% Complementarity index
C = zeros(np,np);
ref = rhst{1};
c = length(ref(~isnan(ref)));
for x = 1:np
    for y = x:np
        C(x,y) = czcmpl(ms3{pci(x)},ms3{pci(y)},c);
        C(y,x) = C(x,y);
    end
end
save compl_index C

cd(mm)

% -------------------------------------------------------------------------
function [irhst,rhst2,thst] = czplace2(ncu,x_pos,y_pos,x_pos_ts)
%CZPLACE2   Place rate map.
%   [IRM RM OM] = CZPLACE2(VD,XPOS,YPOS,TS) returns place rate map (RM), 
%   interpolated place rate map (IRM, for visualization) and occupancy map
%   (OM) for discriminated unit (VD). Input arguments XPOS and YPOS are the
%   rat position coordinates, TS should contain the corresponding timestamps
%   for the position values. Bins with less than 0.25 time or 3 visits are
%   discarded.
%
%   In CZPLACE2, both the occupancy map and the rate map are smoothed by a
%   Gaussian kernel.
%
%   See also CZPLACE.

% Input argumnet check
error(nargchk(4,4,nargin))

% Position vector correction
dx = diff(x_pos);       % leave out the outliers
fdx = find(abs(dx)>8&abs([dx(2:end) 0])>8&abs(dx-[dx(2:end) 0])>16);
x_pos(fdx+1) = [];
y_pos(fdx+1) = [];
x_pos_ts(fdx+1) = [];
dy = diff(y_pos);
fdy = find(abs(dy)>8&abs([dy(2:end) 0])>8&abs(dy-[dy(2:end) 0])>16);
x_pos(fdy+1) = [];
y_pos(fdy+1) = [];
x_pos_ts(fdy+1) = [];

inxs = x_pos > 0 & y_pos > 0;   % leave out (0;0) points
x_pos2 = x_pos(inxs);
y_pos2 = y_pos(inxs);
pos_ts = x_pos_ts(inxs);

% Determine space bins
fminx = min(x_pos2);
fmaxx = max(x_pos2);
fminy = min(y_pos2);
fmaxy = max(y_pos2);

% cst = (fmaxx - fminx) / 25.5;
cst = 8;
xedge = fminx-0.0001:cst:fmaxx+cst;      % spatial bins
yedge = fminy-0.0001:cst:fmaxy+cst;
xbins = length(xedge) - 1;
ybins = length(yedge) - 1;

% Determine the position of action potentials
ncu = ncu(ncu>pos_ts(1)&ncu<pos_ts(end));     % leave action potentials before the 1st pos. timestamp
len = length(ncu);
ap_xpos = zeros(1,len);
ap_ypos = zeros(1,len);
for k = 1:len    % interpolate action potential position
    tap = ncu(k);
    ind = find(pos_ts<tap,1,'last');
    tlow = pos_ts(ind);
    xlow = x_pos2(ind);
    ylow = y_pos2(ind);
    thigh = pos_ts(ind+1);
    xhigh = x_pos2(ind+1);
    yhigh = y_pos2(ind+1);
    mpl = (tap - tlow) / (thigh - tlow);
    ap_xpos(k) = xlow + (xhigh - xlow) * mpl;
    ap_ypos(k) = ylow + (yhigh - ylow) * mpl;
end

hst = hist3([ap_xpos' ap_ypos'],'Edges',{xedge yedge})';     % spike number in space
hst = hst(1:end-1,1:end-1);     % last bin corresponds to the upper edge value

% Calculate the time spent in each space bin
thst = zeros(size(hst));        % time spent in each bin (occupancy map)
tvhst = zeros(size(hst));       % number of visits
for xb = 1:xbins
    for yb = 1:ybins
        xedlow = xedge(xb);
        xedhigh = xedge(xb+1);
        xlin = valuecrossing(pos_ts,x_pos2,xedlow,'up');
        xhin = valuecrossing(pos_ts,x_pos2,xedhigh,'down');
        xlout = valuecrossing(pos_ts,x_pos2,xedlow,'down');
        xhout = valuecrossing(pos_ts,x_pos2,xedhigh,'up');
        xin = union(xlin,xhin);
        xout = union(xlout,xhout);
        if xout(1) < xin(1)
            xin = [pos_ts(1) xin];
        end
        if xin(end) > xout(end)
            xout = [xout pos_ts(end)];
        end
        if ~isequal(length(xin),length(xout))
            error('Technical error 119.')
        end
                
        yedlow = yedge(yb);
        yedhigh = yedge(yb+1);
        ylin = valuecrossing(pos_ts,y_pos2,yedlow,'up');
        yhin = valuecrossing(pos_ts,y_pos2,yedhigh,'down');
        ylout = valuecrossing(pos_ts,y_pos2,yedlow,'down');
        yhout = valuecrossing(pos_ts,y_pos2,yedhigh,'up');
        yin = union(ylin,yhin);
        yout = union(ylout,yhout);
        if yout(1) < yin(1)
            yin = [pos_ts(1) yin];
        end
        if yin(end) > yout(end)
            yout = [yout pos_ts(end)];
        end
        if ~isequal(length(yin),length(yout))
            error('Technical error 137.')
        end
        
        for k1 = 1:length(xin)
            cxin = xin(k1);
            cxout = xout(find(xout>cxin,1,'first'));
            logix = yin<cxout&yout>cxin;
            cyins = yin(logix);
            cyouts = yout(logix);
            for k2 = 1:length(cyins)
                cyin = cyins(k2);
                cyout = cyouts(k2);
                thst(yb,xb) = thst(yb,xb) + (min(cyout,cxout) - max(cyin,cxin));
                tvhst(yb,xb) = tvhst(yb,xb) + 1;
            end
        end
    end
end
gk = gausskernel([3 3]);
thst = smooth2_nonnan(thst,gk); % smoothing
rhst = hst ./ thst;

rhst2 = rhst .* zero2nan(double(thst>=0.25)) .* zero2nan(double(tvhst>=3));       % minimum 3 visits, 0.25 s spent in bin
rhst2 = smooth2_nonnan(rhst2,gk);
irhst = interp2_nonnan(rhst2,5);

% Plot result
% figure
% pcolor(irhst)
% shading flat