function czplaceanalysis2
%CZPLACEANALYSIS2   Analysis of spatial firing.
%   CZPLACEANALYSIS2 seeks for the sign of inhibition on the crosscorrelograms:
%   it compares right central values (3, leaing the central 2 out) with right 
%   side values (6) using Behrens-Fisher test. The program also visualizes
%   the poits of excitation on the rate maps. Edit code to modify input and
%   output directories!
%
%   See also CZPLACE2 and CZPLACEANALYSIS.

% Directories & filename
global DATAPATH
fname = 'hux026-day08-tr3-base-01'
inpdir = [DATAPATH 'Czurko\discriminated2\pos\' fname '\'];
resdir = [DATAPATH 'Czurko\discriminated2\pos\' fname '\'];
fns = findstr(fname,'_');
titlestr = fname;
titlestr(fns) = ' ';
mm = pwd;
cd(resdir)

% Import
load placecell_index
load placemaps
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
placecell = neuron(pci);
np = length(pci);

% Read position data
load([DATAPATH 'Czurko\discriminated2\pos\' fname '\x1_pos.mat']);
x_pos = data';
load([DATAPATH 'Czurko\discriminated2\pos\' fname '\x1_pos_ts.mat']);
x_pos_ts = data';
load([DATAPATH 'Czurko\discriminated2\pos\' fname '\y1_pos.mat']);
y_pos = data';
load([DATAPATH 'Czurko\discriminated2\pos\' fname '\y1_pos_ts.mat']);
y_pos_ts = data';

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

% Inhibition on cross-correlation
for x = 1:np
    for y = x+1:np
        ccr = lxcorr(placecell{x},placecell{y});   % window: +-5 ms
        cnt = (length(ccr) + 1) / 2;
        after1 = ccr(cnt+2:cnt+4);
        after2 = ccr(cnt+5:cnt+10);
        [h,p] = ttest2(after1,after2,0.05,'left','unequal');     % Behrens-Fisher test
        save(['inhib_' num2str(x) '_' num2str(y)],'h','p')
        
% Localization of excitation
        fminx = min(x_pos2);
        fmaxx = max(x_pos2);
        fminy = min(y_pos2);
        fmaxy = max(y_pos2);
        H1 = figure;
        irhst = interp2(rhst{pci(y)},5);
        pirhst = nanpad(irhst,1);
        px = linspace(fminx,fmaxx,size(pirhst,2));
        py = linspace(fminy,fmaxy,size(pirhst,2));
        pcolor(px,py,pirhst)
        A1 = gca;
        shading flat
        hold on
        colormap(bone)
        H2 = figure;
        irhst = interp2(rhst{pci(x)},5);
        pirhst = nanpad(irhst,1);
        pcolor(px,py,pirhst)
        A2 = gca;
        colormap(bone)
        hold on
        shading flat
%         figure;plot(x_pos2,y_pos2);hold on
        for k = 1:length(placecell{y})
            vd_pyr = placecell{y};
            vd_int = placecell{x};
            lvd_int = vd_int(vd_int>vd_pyr(k)&vd_int<vd_pyr(k)+2/1000);
            if ~isempty(lvd_int)
                [ap_xpos ap_ypos] = appos(vd_pyr(k),pos_ts,x_pos2,y_pos2);
%                 figure(H1)
                plot(A1,ap_xpos,ap_ypos,'.','Color','red','MarkerSize',18)
%                 figure(H2)
                plot(A2,ap_xpos,ap_ypos,'.','Color','red','MarkerSize',18)
            end
        end
        title(A1,[num2str(x) '-' num2str(y)])
        title(A2,[num2str(x) '-' num2str(y)])
        str = ['saveas(H1,''locexcit_pyr_' num2str(x) '_' num2str(y) '.fig'')'];
        eval(str)
        str = ['saveas(H2,''locexcit_int_' num2str(x) '_' num2str(y) '.fig'')'];
        eval(str)
        
    end
end

% Complementarity index
% C = zeros(np,np);
% ref = rhst{1};
% c = length(ref(~isnan(ref)));
% for x = 1:np
%     for y = x:np
%         C(x,y) = czcmpl(ms3{pci(x)},ms3{pci(y)},c);
%         C(y,x) = C(x,y);
%     end
% end
% save compl_index C

cd(mm)

% -------------------------------------------------------------------------
function ccr = lxcorr(ncc1,ncc2)
%LXCORR   Crosscorrelation.
%   LXCORR(VD1,VD2) calculates crosscorrelogram for discriminated units
%   VD1 and VD2, using a +-5 ms time window.
%
%   CCR = CZXCORR(VD1,VD2) returns the crosscorrelogram variable.
%
%   See also XCORR and CZACORR.

% Input argument check
error(nargchk(2,2,nargin))

% Calculate spike times in milliseconds
sr = 2000;
nc1 = ncc1 * sr;
nc2 = ncc2 * sr;

% Crosscorrelogram
zunit1 = zeros(1,round(max([nc1; nc2]))+5);
zunit2 = zunit1;
zunit1(round(nc1)) = 1;
zunit2(round(nc2)) = 1;
ccr = xcorr(zunit2,zunit1,0.005*sr);     % 1->2; window: -5 ms - 5 ms
ccr = ccr / length(nc2);     % norm. with the no. of ref. events to get transmission prob.

% -------------------------------------------------------------------------
function [ap_xpos ap_ypos] = appos(tap,pos_ts,x_pos2,y_pos2)

ind = find(pos_ts<tap,1,'last');
tlow = pos_ts(ind);
xlow = x_pos2(ind);
ylow = y_pos2(ind);
thigh = pos_ts(ind+1);
xhigh = x_pos2(ind+1);
yhigh = y_pos2(ind+1);
mpl = (tap - tlow) / (thigh - tlow);
ap_xpos = xlow + (xhigh - xlow) * mpl;
ap_ypos = ylow + (yhigh - ylow) * mpl;