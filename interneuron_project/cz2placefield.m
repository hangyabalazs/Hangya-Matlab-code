function cz2placefield
%CZ2PLACEFIELD   Firing rate map and spatial information.
%   CZ2PLACEFIELD generates smoothed rate maps (5 x 5 pixel Gaussian
%   kernel, sigma = 1; see CZPLACEFIGS2). It calculates and saves spatial
%   information measures: spatial information for all spatial bins, its
%   averaged version and spatial information per spike. See the reference
%   below for details.
%
%   Ref.: Olypher AV, Lansky P, Muller RU, Fenton AA (2003) Quantifying
%   location-specific information in the discharge of rat hippocampal place
%   cells. J Neurosci Meth 127:123-135;
%
%   See also CZPLACEFIGS2.

% Directories
global DATAPATH
inpdir_unit = [DATAPATH 'Czurko2\units2\'];
inpdir_pos = [DATAPATH 'Czurko2\Positions\'];
resdir = [DATAPATH 'Czurko2\Placefields_log2\'];
files = b_filelist(inpdir_unit);
sf = length(files);
mm = pwd;
cd(resdir)

% Progress indicator
wb = waitbar(0,'Running CZ2PLACEFIELD...','Position',[360 250 315 50]);    %Progress indicator
global WB
WB(end+1) = wb;

% Main: calculate place fields
MIPosX = zeros(1,sf);
IPos = zeros(1,sf);
ISpike = zeros(1,sf);
for o = 1:sf
    disp(o)
    fname = files(o).name
    fs = findstr(fname,'_');
    fn = fname(1:fs(1)-1);
    ff = [inpdir_unit fname];    % load unit
    [unit0 unit] = xlsread(ff);
    if isempty(str2num(unit{1}))
        unit(1,:)=[];
        unit = unit(:);
        iu = cellfun(@isempty,unit);
        unit(iu) = [];
    end
    unit = cellfun(@str2num,unit);
    
    [xl1,xl2] = xlsread([inpdir_pos fn '.xls']);    % load position data
    pos = cellfun(@str2num,xl2);
    x_pos_ts = pos(:,1)';
    x_pos_ts = x_pos_ts(~isnan(x_pos_ts));
    x_pos = pos(:,2)';
    x_pos = x_pos(~isnan(x_pos));
    y_pos_ts = pos(:,3)';
    y_pos_ts = y_pos_ts(~isnan(y_pos_ts));
    y_pos = pos(:,4)';
    y_pos = y_pos(~isnan(y_pos));
    
    [irhst_int rhst_int thst Iposx Ipos Ispike] = lczplace(unit,x_pos,y_pos,x_pos_ts); % rate map and spatial inf.
    maxIposx = max(Iposx(:));
    MIPosX(o) = maxIposx;
    IPos(o) = Ipos;
    ISpike(o) = Ispike;
    Hint = figure;     % plot interneuron rate map
    pcolor(nanpad(irhst_int,1))
    shading flat;
    axis off
    colorbar
    fns = [resdir fname(1:end-4) '_RATEMAP.fig'];
    title(['Max. pos. inf.: ' num2str(maxIposx) '   Average pos. inf.: ' num2str(Ipos) ...
        '   Inf. per spike: ' num2str(Ispike)])
    saveas(Hint,fns)    % save
    close(Hint)
    
% Save
    fns = [resdir fname(1:end-4) '_SPINF.mat'];
    save(fns,'Iposx','maxIposx','Ipos','Ispike')
    fns = [resdir fname(1:end-4) '_MAP.mat'];
    save(fns,'irhst_int','rhst_int','thst')
end
fns = [resdir 'spatial_information.mat'];
save(fns,'MIPosX','IPos','ISpike')
close(wb)
cd(mm)

% -------------------------------------------------------------------------
function [irhst,rhst2,thst,Iposx,Ipos,Ispike] = lczplace(ncu,x_pos,y_pos,x_pos_ts)
%LCZPLACE   Place rate map.
%   [IRM RM OM] = LCZPLACE(VD,XPOS,YPOS,TS) returns place rate map (RM), 
%   interpolated place rate map (IRM, for visualization) and occupancy map
%   (OM) for discriminated unit (VD). Input arguments XPOS and YPOS are the
%   rat position coordinates, TS should contain the corresponding timestamps
%   for the position values. Bins with less than 0.25 time or 3 visits are
%   discarded.
%
%   In LCZPLACE, both the occupancy map and the rate map are smoothed by a
%   Gaussian kernel.
%
%   [IRM RM OM IPOSX IPOS ISPIKE] = LCZPLACE(VD,XPOS,YPOS,TS) returns
%   spatial information for all spatial bins, its averaged version and
%   spatial information per spike. See the reference below for details.
%
%   Ref.: Olypher AV, Lansky P, Muller RU, Fenton AA (2003) Quantifying
%   location-specific information in the discharge of rat hippocampal place
%   cells. J Neurosci Meth 127:123-135;
%
%   See also CZPLACE2.

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
figure;plot(x_pos2,y_pos2,'.')

% Determine space bins
fminx = min(x_pos2);
fmaxx = max(x_pos2);
fminy = min(y_pos2);
fmaxy = max(y_pos2);

xedge = fminx-0.0001:8:fmaxx+8;      % spatial bins
yedge = fminy-0.0001:8:fmaxy+8;
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
gk = gausskernel([5 5],1);
thst = smooth2_nonnan(thst,gk); % smoothing
rhst = hst ./ thst;

rhst2 = rhst .* zero2nan(double(thst>=0.25)) .* zero2nan(double(tvhst>=3));       % minimum 3 visits, 0.25 s spent in bin
rhst2 = smooth2_nonnan(rhst2,gk);
irhst = interp2_nonnan(rhst2,5);

% Spatial information
x_center = (xedge(1:end-1) + xedge(2:end)) / 2;
y_center = (yedge(1:end-1) + yedge(2:end)) / 2;
winlen = 0.1;   % window size: 100 ms
maxi = floor((pos_ts(end)-pos_ts(1))/winlen);
K = zeros(1,maxi);
X = zeros(2,maxi);
naninxk = [];
for t = 1:maxi
    inx1 = pos_ts(1) + (t - 1) * winlen;  % Note: overlaping windows!
    inx2 = inx1 + winlen;
    
    lvd = ncu(ncu>=inx1&ncu<inx2);
    K(t) = length(lvd);
    lposx = x_pos2(pos_ts>=inx1&pos_ts<inx2);
    lposy = y_pos2(pos_ts>=inx1&pos_ts<inx2);
    mlposx = mean(lposx);
    mlposy = mean(lposy);
    if isnan(mlposx)
        X(:,t) = [NaN NaN];
        naninxk = [naninxk t];
    else
        X(:,t) = [x_center(abs(x_center-mlposx)==min(abs(x_center-mlposx)))...
            y_center(abs(y_center-mlposy)==min(abs(y_center-mlposy)))];
    end
end
X(isnan(X)) = [];
X = reshape(X,2,length(X)/2);
K(naninxk) = [];
maxk = max(K);
Pk = zeros(1,maxk+1);
for j = 1:maxk+1
    Pk(j) = length(find(K==(j-1))) / length(K);
end
Px = zeros(length(x_center),length(y_center));
for j1 = 1:length(x_center)
    for j2 = 1:length(y_center)
        Px(j1,j2) = length(find(X(1,:)==x_center(j1)&X(2,:)==y_center(j2))) / size(X,2);
    end
end
Pkcx = zeros(maxk+1,length(x_center),length(y_center));
for j = 1:maxk+1
    for j1 = 1:length(x_center)
        for j2 = 1:length(y_center)
            Pkcx(j,j1,j2) = length(find(K==(j-1)&X(1,:)==x_center(j1)&X(2,:)==y_center(j2))) ...
                / length(find(X(1,:)==x_center(j1)&X(2,:)==y_center(j2)));
        end
    end
end
Pkcx = nan2zero(Pkcx);

Iposx = zeros(length(x_center),length(y_center));
for j1 = 1:length(x_center)
    for j2 = 1:length(y_center)
        for j = 1:maxk+1
            ad = Pkcx(j,j1,j2) * log2(Pkcx(j,j1,j2)/Pk(j));
            ad = nan2zero(ad);
            Iposx(j1,j2) = Iposx(j1,j2) + ad;
        end
    end
end

Ipos = 0;
for j1 = 1:length(x_center)
    for j2 = 1:length(y_center)
        Ipos = Ipos + Px(j1,j2) * Iposx(j1,j2);
    end
end

Li = zeros(length(x_center),length(y_center));
for j1 = 1:length(x_center)
    for j2 = 1:length(y_center)
        for j = 1:maxk+1
            Li(j1,j2) = Li(j1,j2) + (j - 1) * Pkcx(j,j1,j2);
        end
    end
end

L = 0;
for j1 = 1:length(x_center)
    for j2 = 1:length(y_center)
        ad = Li(j1,j2) * Px(j1,j2);
        ad = nan2zero(ad);
        L = L + ad;
    end
end

Ispike = 0;
for j1 = 1:length(x_center)
    for j2 = 1:length(y_center)
        ad = Li(j1,j2) * log2(Li(j1,j2)/L) * Px(j1,j2);
        ad = nan2zero(ad);
        Ispike = Ispike + ad;
    end
end
Ispike = (1 / L) * Ispike;

% Plot result
% figure
% pcolor(irhst)
% shading flat