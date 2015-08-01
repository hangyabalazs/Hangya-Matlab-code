function [irhst, rhst, rhst2, thst, thst2] = ...
    place(spike_times,x_pos,y_pos,x_pos_ts,varargin)
%PLACE   Spatial firing map.
%   [IRM RM RM2 OM OM2] = PLACE(VD,XPOS,YPOS,TS) calculates spatial firing
%   maps. Spatial bins with less than 0.25s cumulative dwell time or less 
%   than 3 visits are discarded. Both the occupancy map and the rate map 
%   are smoothed by a Gaussian kernel over 3x3 bins.
%
%   Input arguments:
%       VD - spike times (in seconds)
%       XPOS - X position coordinates
%       YPOS - Y position coordinates
%       TS - time stamps correponding to the position coordinates
%
%   Output arguments:
%       IRM - interpolated rate map (for visualization)
%       RM - raw rate map
%       RM2 - rate map smoothed by a Gaussian kernel
%       OM - raw occupancy map
%       OM2 - occupancy map smoothed by a Gaussian kernel
%
%   Optional input parameter-value pairs, with default values:
%       BinSize, 8 - size of the spatial bins
%       Display, false - controls plotting of the rate map

%   Balazs Hangya
%   Cold Spring Harbor Laboratory
%   1 Bungtown Road, Cold Spring Harbor
%   New York 11724, USA
%   balazs.cshl@gmail.com
%
%   Institute of Experimental Medicine
%   Szigony street 43, Budapest 1083, Hungary
%   hangyab@koki.hu

% Input argumnet check
prs = inputParser;
addRequired(prs,'spike_times',@isnumeric)   % spike times
addRequired(prs,'x_pos',@isnumeric)   % X coordinates
addRequired(prs,'y_pos',@isnumeric)   % Y coordinates
addRequired(prs,'x_pos_ts',@isnumeric)  % time stamps for the coordinates
addParamValue(prs,'BinSize',8,@isnumeric)   % spatial bin size (default: 8)
addParamValue(prs,'Display',false,@(s)islogical(s)|ismember(s,[0 1]))   % control displaying
parse(prs,spike_times,x_pos,y_pos,x_pos_ts,varargin{:})
g = prs.Results;

% Position vector correction
x_pos = x_pos(:)';
y_pos = y_pos(:)';
x_pos_ts = x_pos_ts(:)';
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
fminx = min(x_pos2);   % boundaries
fmaxx = max(x_pos2);
fminy = min(y_pos2);
fmaxy = max(y_pos2);

xedge = fminx-0.0001:g.BinSize:fmaxx+g.BinSize;      % spatial bins
yedge = fminy-0.0001:g.BinSize:fmaxy+g.BinSize;
xbins = length(xedge) - 1;
ybins = length(yedge) - 1;

% Determine the position of action potentials
spike_times = spike_times(spike_times>pos_ts(1)&spike_times<pos_ts(end));     % leave action potentials before the 1st pos. timestamp
len = length(spike_times);   % number of spikes
ap_xpos = zeros(1,len);
ap_ypos = zeros(1,len);
for k = 1:len    % interpolate action potential position
    tap = spike_times(k);   % spike time
    ind = find(pos_ts<tap,1,'last');
    tlow = pos_ts(ind);   % corresponding position time stamp (before)
    xlow = x_pos2(ind);   % corresponding x coordinate (before)
    ylow = y_pos2(ind);   % corresponding y coordinate (before)
    thigh = pos_ts(ind+1);   % corresponding position time stamp (after)
    xhigh = x_pos2(ind+1);   % corresponding x coordinate (after)
    yhigh = y_pos2(ind+1);   % corresponding y coordinate (after)
    mpl = (tap - tlow) / (thigh - tlow);   % interpolate
    ap_xpos(k) = xlow + (xhigh - xlow) * mpl;
    ap_ypos(k) = ylow + (yhigh - ylow) * mpl;
end

hst = hist3([ap_xpos' ap_ypos'],'Edges',{xedge yedge})';     % spike number in space
hst = hst(1:end-1,1:end-1);     % last bin corresponds to the upper edge value

% Calculate the time spent in each space bin (occupancy)
thst = zeros(size(hst));        % time spent in each bin (occupancy map)
tvhst = zeros(size(hst));       % number of visits
for xb = 1:xbins   % loop through spatial bins
    for yb = 1:ybins
        xedlow = xedge(xb);
        xedhigh = xedge(xb+1);
        xlin = valuecrossing(pos_ts,x_pos2,xedlow,'up');   % x entry
        xhin = valuecrossing(pos_ts,x_pos2,xedhigh,'down');
        xlout = valuecrossing(pos_ts,x_pos2,xedlow,'down');   % x exit
        xhout = valuecrossing(pos_ts,x_pos2,xedhigh,'up');
        xin = union(xlin,xhin);  % x entry
        xout = union(xlout,xhout);   % x exit
        if (isempty(xin) && ~isempty(xout)) || xout(1) < xin(1)
            xin = [pos_ts(1) xin];   % first bin
        end
        if (~isempty(xin) && isempty(xout)) || xin(end) > xout(end)
            xout = [xout pos_ts(end)];   % last bin
        end
        if ~isequal(length(xin),length(xout))
            error('Technical error 119.')
        end
                
        yedlow = yedge(yb);
        yedhigh = yedge(yb+1);
        ylin = valuecrossing(pos_ts,y_pos2,yedlow,'up');   % y entry
        yhin = valuecrossing(pos_ts,y_pos2,yedhigh,'down');
        ylout = valuecrossing(pos_ts,y_pos2,yedlow,'down');   % y exit
        yhout = valuecrossing(pos_ts,y_pos2,yedhigh,'up');
        yin = union(ylin,yhin);   % y entry
        yout = union(ylout,yhout);   % y exit
        if (isempty(yin) && ~isempty(yout)) || yout(1) < yin(1)
            yin = [pos_ts(1) yin];   % first bin
        end
        if (~isempty(yin) && isempty(yout)) || yin(end) > yout(end)
            yout = [yout pos_ts(end)];   % last bin
        end
        if ~isequal(length(yin),length(yout))
            error('Technical error 137.')
        end
        
        for k1 = 1:length(xin)   % construct occupancy map
            cxin = xin(k1);
            cxout = xout(find(xout>cxin,1,'first'));
            logix = yin<cxout&yout>cxin;
            cyins = yin(logix);
            cyouts = yout(logix);
            for k2 = 1:length(cyins)
                cyin = cyins(k2);
                cyout = cyouts(k2);
                thst(yb,xb) = thst(yb,xb) + (min(cyout,cxout) - max(cyin,cxin));   % time in bin
                tvhst(yb,xb) = tvhst(yb,xb) + 1;   % number of visits
            end
        end
    end
end

% Rate map
gk = gausskernel([3 3]);   % smoothing kernel
thst2 = smooth2_nonnan(thst,gk);   % smoothed occupancy map
rhst = hst ./ thst2;   % raw rate map
rhst2 = rhst .* zero2nan(double(thst2>=0.25)) .* zero2nan(double(tvhst>=3));   % apply minimal occupancy criteria       % minimum 3 visits, 0.25 s spent in bin
rhst2 = smooth2_nonnan(rhst2,gk);   % smoothed rate map
irhst = interp2_nonnan(rhst2,5);   % interpolated rate map

% Plot rate map
if g.Display
    figure
    pcolor(irhst)
    shading flat
end