function czplaceanalysis3b(fname)
%CZPLACEANALYSIS2B   Analysis of spatial firing.
%   CZPLACEANALYSIS2B calculates the temporal distribution of monosynaptic
%   excitation, the distance of the points of excitation from the pyr. cell
%   place field and crosscorrelation inside/outside the pyramidal cell's 
%   place field. Edit code to modify input and output directories!
%
%   See also CZPLACEANALYSIS and CZPLACEANALYSIS2B.

% Directories & filename
global DATAPATH
nop = 'neg';
inpdir = [DATAPATH 'Czurko\discriminated2\' nop '\' fname '\'];
resdir = [DATAPATH 'Czurko\discriminated2\' nop '\' fname '\'];
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
load([DATAPATH 'Czurko\discriminated2\' nop '\' fname '\x1_pos.mat']);
x_pos = data';
load([DATAPATH 'Czurko\discriminated2\' nop '\' fname '\x1_pos_ts.mat']);
x_pos_ts = data';
load([DATAPATH 'Czurko\discriminated2\' nop '\' fname '\y1_pos.mat']);
y_pos = data';
load([DATAPATH 'Czurko\discriminated2\' nop '\' fname '\y1_pos_ts.mat']);
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

% Main
for x = 1:np
    for y = x+1:np
        fminx = min(x_pos2);
        fmaxx = max(x_pos2);
        fminy = min(y_pos2);
        fmaxy = max(y_pos2);
        xedge = fminx-0.0001:8:fmaxx+8;      % spatial bins
        yedge = fminy-0.0001:8:fmaxy+8;
        xbins = length(xedge) - 1;
        ybins = length(yedge) - 1;
        xbinstep = xedge(2) - xedge(1);
        ybinstep = yedge(2) - yedge(1);
        bincenter_xpos = (xedge(1:end-1) + xedge(2:end)) / 2;
        bincenter_ypos = (yedge(1:end-1) + yedge(2:end)) / 2;
        rm_pyr = rhst{pci(y)};
        ri_pyr = rm_pyr .* (double(rm_pyr>b_max_nonnan(rm_pyr)*0.5));
        vd_pyr = placecell{y};
        vd_int = placecell{x};
        vd_pyr = vd_pyr(vd_pyr>pos_ts(1)&vd_pyr<pos_ts(end));     % leave action potentials before the 1st pos.
        vd_int = vd_int(vd_int>pos_ts(1)&vd_int<pos_ts(end));
        vd_trm = [];
        vd_pyr_infield = [];
        vd_pyr_outfield = [];
        dst_excit = [];
        dst_all = [];
        for k = 1:length(vd_pyr)
            lvd_int = vd_int(vd_int>vd_pyr(k)&vd_int<vd_pyr(k)+2/1000);
            [ap_xpos ap_ypos] = appos(vd_pyr(k),pos_ts,x_pos2,y_pos2);
            ap_xbin = floor((ap_xpos-xedge(1))/xbinstep) + 1;
            ap_ybin = floor((ap_ypos-yedge(1))/ybinstep) + 1;
            if ~isempty(lvd_int)
                vd_trm(end+1) = vd_pyr(k);      % time course of excitation
                dst = [];
                for t = 1:numel(ri_pyr)     % distance of excitation from the place field
                    if ~isnan(ri_pyr(t)) && ri_pyr(t) ~= 0
                        [posiy posix] = ind2sub(size(ri_pyr),t);
                        dst(end+1) = sqrt((ap_xpos-bincenter_xpos(posix)).^2+(ap_ypos-bincenter_ypos(posiy)).^2);
                    end
                end
                dst_excit = [dst_excit min(dst)];
                dst_all = [dst_all min(dst)];
            else
                dst = [];
                for t = 1:numel(ri_pyr)
                    if ~isnan(ri_pyr(t)) && ri_pyr(t) ~= 0
                        [posiy posix] = ind2sub(size(ri_pyr),t);
                        dst(end+1) = sqrt((ap_xpos-bincenter_xpos(posix)).^2+(ap_ypos-bincenter_ypos(posiy)).^2);
                    end
                end
                dst_all = [dst_all min(dst)];
            end
            if ~isnan(ri_pyr(ap_ybin,ap_xbin)) && ri_pyr(ap_ybin,ap_xbin) ~= 0
                vd_pyr_infield = [vd_pyr_infield vd_pyr(k)];
            else
                vd_pyr_outfield = [vd_pyr_outfield vd_pyr(k)];
            end
        end
        vd_int_infield = [];
        vd_int_outfield = [];
        for k = 1:length(vd_int)
            [ap_xpos ap_ypos] = appos(vd_int(k),pos_ts,x_pos2,y_pos2);
            ap_xbin = floor((ap_xpos-xedge(1))/xbinstep) + 1;
            ap_ybin = floor((ap_ypos-yedge(1))/ybinstep) + 1;
            if ~isnan(ri_pyr(ap_ybin,ap_xbin)) && ri_pyr(ap_ybin,ap_xbin) ~= 0
                vd_int_infield = [vd_int_infield vd_int(k)];
            else
                vd_int_outfield = [vd_int_outfield vd_int(k)];
            end
        end
        
        vd_trm_infield = [];
        rn = rand(size(vd_trm)) / 50;
        rn2 = [];
        for k = 1:length(vd_pyr_infield)
            lvd_int = vd_int(vd_int>vd_pyr_infield(k)&vd_int<vd_pyr_infield(k)+2/1000);
            if ~isempty(lvd_int)
                vd_trm_infield(end+1) = vd_pyr_infield(k);      % time of infield excitation
                fn = find(vd_trm==vd_pyr_infield(k),1,'first');
                rn2(end+1) = rn(fn);
            end
        end
        if ~isempty(vd_trm)
            vd_pyr1000 = round(vd_pyr*1000);
            vd_trm1000 = round(vd_trm*1000);
            vd_trm_infield1000 = round(vd_trm_infield*1000);
            z = zeros(1,vd_pyr1000(end));
            z(vd_trm1000) = 1 + rn;
            z_infield = zeros(1,vd_pyr1000(end));
            z_infield(vd_trm_infield1000) = 1 + rn2;
            figure
            plot(linspace(0,vd_trm1000(end),length(z))/1000,z,'.')
            ylim([0 2])
            hold on
            if ~isempty(vd_trm_infield)
                plot(linspace(0,vd_trm1000(end),length(z_infield))/1000,z_infield,'r.')
            end
            title([num2str(x) '-' num2str(y)])
            str = ['saveas(gcf,''timeexcit_' num2str(x) '_' num2str(y) '.fig'')'];
            eval(str)   % save (time of excitation)
            zt = z(1:floor(length(z)/100)*100);
            z2 = reshape(zt,length(zt)/100,100);
            figure
            bar(1:100,sum(z2))
            title([num2str(x) '-' num2str(y)])
            str = ['saveas(gcf,''timeexcit2_' num2str(x) '_' num2str(y) '.fig'')'];
            eval(str)
            close all
            mdst_all = mean(dst_all);   % save (distance of excitation)
            mdst_excit = mean(dst_excit);
            str = ['save(''excitdist_' num2str(x) '_' num2str(y)...
                ''',''dst_all'',''dst_excit'',''mdst_all'',''mdst_excit'')'];
            eval(str)
        end
        
        % infield and outfield crosscorrelation
%         [H1 H2 trsc_infield] = czxcorr(vd_int_infield',vd_pyr_infield');   % window: +-50 ms
%         title([num2str(x) '-' num2str(y)])
%         str = ['saveas(gcf,''normcrosscorr_infield_' num2str(x) '_' num2str(y) '.fig'')'];  % save
%         eval(str)
%         str = ['saveas(gcf,''normcrosscorr_infield_' num2str(x) '_' num2str(y) '.jpg'')'];
%         eval(str)
%         figure(H1)
%         title([num2str(x) '-' num2str(y)])
%         str = ['saveas(gcf,''crosscorr_infield_' num2str(x) '_' num2str(y) '.fig'')'];
%         eval(str)
%         str = ['saveas(gcf,''crosscorr_infield_' num2str(x) '_' num2str(y) '.jpg'')'];
%         eval(str)
%         str = ['save(''trsc_infield_' num2str(x) '_' num2str(y) ''',''trsc_infield'')'];
%         eval(str)
%         close all
%         
%         [H1 H2] = czxcorr2(vd_int_infield',vd_pyr_infield');   % window: +-5 ms
%         title([num2str(x) '-' num2str(y)])
%         str = ['saveas(gcf,''nsmallxcorr_infield_' num2str(x) '_' num2str(y) '.fig'')'];  % save
%         eval(str)
%         str = ['saveas(gcf,''nsmallxcorr_infield_' num2str(x) '_' num2str(y) '.jpg'')'];
%         eval(str)
%         figure(H1)
%         title([num2str(x) '-' num2str(y)])
%         str = ['saveas(gcf,''smallxcorr_infield' num2str(x) '_' num2str(y) '.fig'')'];
%         eval(str)
%         str = ['saveas(gcf,''smallxcorr_infield' num2str(x) '_' num2str(y) '.jpg'')'];
%         eval(str)
%         close all
%         
%         [H1 H2 trsc_outfield] = czxcorr(vd_int_outfield',vd_pyr_outfield');   % window: +-50 ms
%         title([num2str(x) '-' num2str(y)])
%         str = ['saveas(gcf,''normcrosscorr_outfield_' num2str(x) '_' num2str(y) '.fig'')'];  % save
%         eval(str)
%         str = ['saveas(gcf,''normcrosscorr_outfield_' num2str(x) '_' num2str(y) '.jpg'')'];
%         eval(str)
%         figure(H1)
%         title([num2str(x) '-' num2str(y)])
%         str = ['saveas(gcf,''crosscorr_outfield_' num2str(x) '_' num2str(y) '.fig'')'];
%         eval(str)
%         str = ['saveas(gcf,''crosscorr_outfield_' num2str(x) '_' num2str(y) '.jpg'')'];
%         eval(str)
%         str = ['save(''trsc_outfield_' num2str(x) '_' num2str(y) ''',''trsc_outfield'')'];
%         eval(str)
%         close all
%         
%         [H1 H2] = czxcorr2(vd_int_outfield',vd_pyr_outfield');   % window: +-5 ms
%         title([num2str(x) '-' num2str(y)])
%         str = ['saveas(gcf,''nsmallxcorr_outfield_' num2str(x) '_' num2str(y) '.fig'')'];  % save
%         eval(str)
%         str = ['saveas(gcf,''nsmallxcorr_outfield_' num2str(x) '_' num2str(y) '.jpg'')'];
%         eval(str)
%         figure(H1)
%         title([num2str(x) '-' num2str(y)])
%         str = ['saveas(gcf,''smallxcorr_outfield' num2str(x) '_' num2str(y) '.fig'')'];
%         eval(str)
%         str = ['saveas(gcf,''smallxcorr_outfield' num2str(x) '_' num2str(y) '.jpg'')'];
%         eval(str)
        close all
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