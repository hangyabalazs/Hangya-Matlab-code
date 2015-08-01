function [C Cb Cc trsr] = czplacefigs3(int,pyr,exp,type)
%CZPLACEFIGS3   Prepare figures for publication.
%   CZPLACEFIGS3 generates smoothed rate maps (5 x 5 pixel Gaussian kernel,
%   sigma = 1), black auto- and crosscorrelograms; new complementarity
%   indeces (see CZCMPL2) and transmission success rate values (sum of
%   three bins of normalized crosscorrelogram from 1 to 3 ms).
%
%   Syntax:
%   [C CB CC TRSR] = CZPLACEFIGS3(INT,PYR,EXP,TYPE)
%
%   Input arguments:
%   INT: place cell index of the interneuron,
%   PYR: place cell index of the pyramidal cell,
%   EXP: name of the animals directory,
%   TYPE: type of the pair ('neg') or ('pos').
%
%   Output arguments:
%   C: type I complementarity index,
%   CB: type II complementarity index,
%   CC: type III complementarity index,
%   TRSR: transmission success rate.
%
%   See also CZPLACEANALYSIS and CZCMPL2.

% Load previous data
global DATAPATH
% type = 'new'
inpdir = [DATAPATH 'Czurko\discriminated2\' type '\' exp '\'];
if ~isdir(inpdir)
    inpdir = [DATAPATH 'Czurko\discriminated2\new\' exp '\'];
end
resdir = [inpdir 'final2\'];
if ~isdir(resdir)
    mkdir(resdir)
end
cd(inpdir)
load('pfs.mat')
load('placemaps.mat')
load('compl_index.mat')
load('placecell_index.mat')
int2 = find(pci==int);
pyr2 = find(pci==pyr);

% Import unit
fl = dir(inpdir);   % find text files
txts = {};
for k = 3:length(fl)
    if isequal(fl(k).name(end-2:end),'mat') && isequal(fl(k).name(1:2),'nr')
        txts{end+1} = fl(k).name;
    end
end
load([inpdir txts{int}])
neuron_int = data;
load([inpdir txts{pyr}])
neuron_pyr = data;

% Read position data
load([inpdir '\x1_pos.mat']);
x_pos = data';
load([inpdir '\x1_pos_ts.mat']);
x_pos_ts = data';
load([inpdir '\y1_pos.mat']);
y_pos = data';
load([inpdir '\y1_pos_ts.mat']);
y_pos_ts = data';

% Rate maps
[irhst_int rhst_int thst] = lczplace(neuron_int,x_pos,y_pos,x_pos_ts);
[irhst_pyr rhst_pyr thst] = lczplace(neuron_pyr,x_pos,y_pos,x_pos_ts);
Hint = figure;     % plot interneuron rate map
pcolor(nanpad(irhst_int,1))
shading flat;
grid on
colorbar
set(gca,'Layer','top','LineWidth',2,'XColor',[0.3 0.3 0.3],'YColor',...
    [0.3 0.3 0.3],'XTickLabel',{},'YTickLabel',{})
fns = [resdir 'RATEMAP_int' num2str(int) '.fig'];
saveas(Hint,fns)    % save
Hpyr = figure;      % plot pyramidal cell rate map
pcolor(nanpad(irhst_pyr,1))
shading flat;
grid on
colorbar
set(gca,'Layer','top','LineWidth',2,'XColor',[0.3 0.3 0.3],'YColor',...
    [0.3 0.3 0.3],'XTickLabel',{},'YTickLabel',{})
fns = [resdir 'RATEMAP_pyr' num2str(pyr) '.fig'];
saveas(Hpyr,fns)    % save

% mn_int = (min(irhst_int(:)));
% mx_int = (max(irhst_int(:)));
% mn_pyr = (min(irhst_pyr(:)));
% mx_pyr = (max(irhst_pyr(:)));
% figure
% contour(irhst_int,mn_int+0.33*(mx_int-mn_int))
% hold on
% contour(irhst_pyr,mn_pyr+0.67*(mx_pyr-mn_pyr))
% 
% irhst_intn = irhst_int / nansum(irhst_int(:));
% irhst_pyrn = irhst_pyr / nansum(irhst_pyr(:));
% figure
% pcolor(nanpad(irhst_pyrn+irhst_intn,1));shading flat;

% Complementarity index
s_int = rhst_int;
ms3_int = s_int .* (zero2nan(double(s_int>b_max_nonnan(s_int)*0.5)));
s_pyr = rhst_pyr;
ms3_pyr = s_pyr .* (zero2nan(double(s_pyr>b_max_nonnan(s_pyr)*0.5)));
c = length(rhst_int(~isnan(rhst_int)));
[C Cb Cc] = czcmpl2(ms3_int,ms3_pyr,c);

% Auto- and crosscorrelation plots, transmission success rate
open(['autocorr' num2str(int2) '.fig'])
A_autoint = gca;
open(['autocorr' num2str(pyr2) '.fig'])
A_autopyr = gca;
open(['crosscorr_' num2str(min(pyr2,int2)) '_' num2str(max(pyr2,int2)) '.fig'])
A_cross = gca;
open(['normcrosscorr_' num2str(min(pyr2,int2)) '_' num2str(max(pyr2,int2)) '.fig'])
A_normcross = gca;
open(['smallxcorr_' num2str(min(pyr2,int2)) '_' num2str(max(pyr2,int2)) '.fig'])
A_smallcross = gca;

ach1 = findobj(allchild(A_autoint),'Type','line');  % interneuron autocorrelation
ach2 = findobj(allchild(A_autoint),'Type','hggroup');
set(ach2,'FaceColor',[0 0 0],'BarWidth',1)
set(A_autoint,'XLim',[-200 200])
axes(A_autoint)
axis off
fns = [resdir 'AUTO_int' num2str(int) '.fig'];
saveas(gcf,fns)    % save

ach1 = findobj(allchild(A_autopyr),'Type','line');  % place cell autocorrelation
ach2 = findobj(allchild(A_autopyr),'Type','hggroup');
set(ach2,'FaceColor',[0 0 0],'BarWidth',1)
set(A_autopyr,'XLim',[-200 200])
axes(A_autopyr)
axis off
fns = [resdir 'AUTO_pyr' num2str(pyr) '.fig'];
saveas(gcf,fns)    % save

ach1 = findobj(allchild(A_cross),'Type','line');    % crosscorrelation
ach2 = findobj(allchild(A_cross),'Type','hggroup');
set(ach2,'FaceColor',[0 0 0],'BarWidth',1)
set(A_cross,'XLim',[-50 50])
xd = get(ach2,'XData');
yd = get(ach2,'YData');
axes(A_cross)
axis off
H_cross = figure;
bar(xd,fliplr(yd),'FaceColor',[0 0 0],'BarWidth',1)
set(gca,'XLim',[-50 50])
axis off
fns = [resdir 'CROSS_int' num2str(int) '_pyr' num2str(pyr) '.fig'];
saveas(gcf,fns)    % save

ach1 = findobj(allchild(A_normcross),'Type','line');    % normalized crosscorrelation
ach2 = findobj(allchild(A_normcross),'Type','hggroup');
set(ach2,'FaceColor',[0 0 0],'BarWidth',1)
yd = get(ach2,'YData');
trsr = yd(48) + yd(49) + yd(50);
set(A_normcross,'XLim',[-50 50])
xd = get(ach2,'XData');
yd = get(ach2,'YData');
xdl = get(ach1(3),'XData');
ydl = get(ach1(3),'YData');
axes(A_normcross)
axis off
H_normcross = figure;
bar(xd,fliplr(yd),'FaceColor',[0 0 0],'BarWidth',1)
line(xdl,ydl,'Color','r','LineWidth',2)
set(gca,'XLim',[-50 50])
axis off
fns = [resdir 'NCROSS_int' num2str(int) '_pyr' num2str(pyr) '.fig'];
saveas(gcf,fns)    % save

ach1 = findobj(allchild(A_smallcross),'Type','line');    % crosscorrelation
ach2 = findobj(allchild(A_smallcross),'Type','hggroup');
set(ach2,'FaceColor',[0 0 0],'BarWidth',1)
set(A_smallcross,'XLim',[-5 5])
xd = get(ach2,'XData');
yd = get(ach2,'YData');
axes(A_smallcross)
axis off
H_smallcross = figure;
bar(xd,fliplr(yd),'FaceColor',[0 0 0],'BarWidth',1)
set(gca,'XLim',[-5 5])
axis off
fns = [resdir 'SCROSS_int' num2str(int) '_pyr' num2str(pyr) '.fig'];
saveas(gcf,fns)    % save
close all



% -------------------------------------------------------------------------
function [irhst,rhst2,thst] = lczplace(ncu,x_pos,y_pos,x_pos_ts)
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
gk = gausskernel([3 3],1);
thst = smooth2_nonnan(thst,gk); % smoothing
rhst = hst ./ thst;

rhst2 = rhst .* zero2nan(double(thst>=0.25)) .* zero2nan(double(tvhst>=3));       % minimum 3 visits, 0.25 s spent in bin
rhst2 = smooth2_nonnan(rhst2,gk);
irhst = interp2_nonnan(rhst2,5);

% Plot result
figure
pcolor(rhst)
sh = findobj(allchild(gcf),'Type','surface');
set(sh,'EdgeColor','white');
c2 = [1 1 0.502; 0.9843 0.8588 0.4745; 1 0 0; 0 0.651 0; 0.502 0.502 1; 0.4157 0 0.8353];
colormap(c2)