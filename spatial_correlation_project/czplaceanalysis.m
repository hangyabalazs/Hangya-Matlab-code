function czplaceanalysis
%CZPLACEANALYSIS   Analysis of spatial firing.
%   CZPLACEANALYSIS calculates place rate maps, autocorrelation,
%   crosscorrelation, place field similarity and complementarity index for
%   neuronal data. Edit code to modify input and output directories!
%
%   See also CZXCORR, CZACORR, CZPLACE2, CZPFSIZE, CZPFS and CZCMPL.

% Directories & filename
global DATAPATH
fname = 'acin05s032g'
inpdir = [DATAPATH 'Czurko\discriminated2\new3\' fname '\'];
resdir = [DATAPATH 'Czurko\discriminated2\new3\' fname '\'];
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
load([DATAPATH 'Czurko\discriminated2\new3\' fname '\x1_pos.mat']);
x_pos = data';
load([DATAPATH 'Czurko\discriminated2\new3\' fname '\x1_pos_ts.mat']);
x_pos_ts = data';
load([DATAPATH 'Czurko\discriminated2\new3\' fname '\y1_pos.mat']);
y_pos = data';
load([DATAPATH 'Czurko\discriminated2\new3\' fname '\y1_pos_ts.mat']);
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