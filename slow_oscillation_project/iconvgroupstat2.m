function iconvgroupstat2
%ICONCVGROUPSTAT2   Group statistics for convergence and divergence maps.
%   ICONVGROUPSTAT2 calculates convergence and divergence strength
%   distributions as well as distance-contribution strength functions.
%   Imported variables are generated by ICONTRIBMAP.
%
%   See also ICONTRIBMAP.

% Read Excel file
global DATAPATH
fn = [DATAPATH 'Ulbert\patients_all.xls'];
headerrows = 0;
[mtx ntx atx] = xlsread(fn);
ntx(1:headerrows,:) = [];
atx(1:headerrows,:) = [];

% Call
conv_strength = [];
div_strength = [];
conv_strength_rand = [];
div_strength_rand = [];
conv_strength_norm = [];
div_strength_norm = [];
conv_strength_rand_norm = [];
div_strength_rand_norm = [];
conv_strength1 = [];
div_strength1 = [];
conv_strength_rand1 = [];
div_strength_rand1 = [];
conv_strength2 = [];
div_strength2 = [];
conv_strength_rand2 = [];
div_strength_rand2 = [];
conv_strength3 = [];
div_strength3 = [];
conv_strength_rand3 = [];
div_strength_rand3 = [];
dist_conv = [];
dist_div = [];
where_from = [];
where_to = [];
T = 1:26;
for k = T
    k
    pat = atx{k,1};
    patno = num2str(atx{k,2});
    eg = num2str(atx{k,3});
    nm_rows = atx{k,4};
    nm_cols = atx{k,5};
    
    [sCnv sDiv convdm divdm wf wt] = lload(pat,patno,eg,nm_rows,nm_cols);
    [sCnv_rand sDiv_rand] = lload2(pat,patno,eg,nm_rows,nm_cols);
    sCnv = sCnv';
    sDiv = sDiv';
    sCnv_rand = sCnv_rand';
    sDiv_rand = sDiv_rand';
    switch mod(k,3)
        case 1
            conv_strength1 = [conv_strength1 sCnv(1:16)];
            div_strength1 = [div_strength1 sDiv(1:16)];
            conv_strength_rand1 = [conv_strength_rand1 sCnv(1:16)];
            div_strength_rand1 = [div_strength_rand1 sDiv(1:16)];
        case 2
            conv_strength2 = [conv_strength2 sCnv(1:16)];
            div_strength2 = [div_strength2 sDiv(1:16)];
            conv_strength_rand2 = [conv_strength_rand2 sCnv(1:16)];
            div_strength_rand2 = [div_strength_rand2 sDiv(1:16)];
        case 0
            conv_strength3 = [conv_strength3 sCnv(1:16)];
            div_strength3 = [div_strength3 sDiv(1:16)];
            conv_strength_rand3 = [conv_strength_rand3 sCnv(1:16)];
            div_strength_rand3 = [div_strength_rand3 sDiv(1:16)];
    end
    conv_strength = [conv_strength sCnv(1:16)];
    div_strength = [div_strength sDiv(1:16)];
    conv_strength_rand = [conv_strength_rand sCnv_rand(1:16)];
    div_strength_rand = [div_strength_rand sDiv_rand(1:16)];
    conv_strength_norm = [conv_strength_norm (sCnv(1:16)/60000-mean(sCnv_rand(1:16)/60000))/std(sCnv_rand(1:16)/60000)];
    div_strength_norm = [div_strength_norm (sDiv(1:16)/1200-mean(sDiv_rand(1:16)/1200))/std(sDiv_rand(1:16)/1200)];
    conv_strength_rand_norm = [conv_strength_rand_norm (sCnv_rand(1:16)/60000-mean(sCnv_rand(1:16)/60000))/std(sCnv_rand(1:16)/60000)];
    div_strength_rand_norm = [div_strength_rand_norm (sDiv_rand(1:16)/1200-mean(sDiv_rand(1:16)/1200))/std(sDiv_rand(1:16)/1200)];
    dist_conv = [dist_conv convdm];
    dist_div = [dist_div divdm];
    where_from = [where_from wf/sum(wf)];
    where_to = [where_to wt/sum(wt)];
    close all
end

I = conv_strength > 0 | div_strength > 0;
conv_strength = conv_strength(I);
div_strength = div_strength(I);
conv_strength_rand = conv_strength_rand(conv_strength_rand>0);
div_strength_rand = div_strength_rand(div_strength_rand>0);
conv_strength1 = conv_strength1(conv_strength1>0);
div_strength1 = div_strength1(div_strength1>0);
conv_strength_rand1 = conv_strength_rand1(conv_strength_rand1>0);
div_strength_rand1 = div_strength_rand1(div_strength_rand1>0);
conv_strength2 = conv_strength2(conv_strength2>0);
div_strength2 = div_strength2(div_strength2>0);
conv_strength_rand2 = conv_strength_rand2(conv_strength_rand2>0);
div_strength_rand2 = div_strength_rand2(div_strength_rand2>0);
conv_strength3 = conv_strength3(conv_strength3>0);
div_strength3 = div_strength3(div_strength3>0);
conv_strength_rand3 = conv_strength_rand3(conv_strength_rand3>0);
div_strength_rand3 = div_strength_rand3(div_strength_rand3>0);
conv_strength_norm = conv_strength_norm(I);
div_strength_norm = div_strength_norm(I);
conv_strength_rand_norm = conv_strength_rand_norm(I);
div_strength_rand_norm = div_strength_rand_norm(I);

% Strength of convergence in distance bins
udc = unique(dist_conv);
udc(1) = [];
bwf = cell(1,size(udc,2));
for n = 1:size(udc,2)
    bwf{n} = where_from(find(dist_conv==udc(n)));
    if size(bwf{n},2)>= 4
%         mwf(n)= median(bwf{n}); % mean per distance bin of where_from
%         p25wf(n) = prctile(bwf{n},25); % 0.25 percentile per distance bin of where_from
%         p75wf(n) = prctile(bwf{n},75); % 0.75 percentile per distance bin of where_from
        mwf(n) = mean(bwf{n});
        p25wf(n) = mean(bwf{n}) - std(bwf{n}) / sqrt(length(bwf{n}));
        p75wf(n) = mean(bwf{n}) + std(bwf{n}) / sqrt(length(bwf{n}));
    end
end
udc = udc(1:size(mwf,2));
figure; errorbar(udc, mwf, mwf-p25wf, p75wf-mwf, ':'); hold on
X = [ones(size(udc')) udc'];
[b,bint,r,rint,stats]= regress(  mwf', X);
yfit = b(1)+ b(2)* udc;
plot(udc, yfit, 'r')
clear bwf mwf p25wf p75wf

% Strength of divergence in distance bins
udc = unique(dist_div);
udc(1) = [];
bwf = cell(1,size(udc,2));
for n = 1:size(udc,2)
    bwf{n} = where_to(find(dist_div==udc(n)));
    if size(bwf{n},2)>= 4
%         mwf(n)= median(bwf{n}); % mean per distance bin of where_from
%         p25wf(n) = prctile(bwf{n},25); % 0.25 percentile per distance bin of where_from
%         p75wf(n) = prctile(bwf{n},75); % 0.75 percentile per distance bin of where_from
        mwf(n) = mean(bwf{n});
        p25wf(n) = mean(bwf{n}) - std(bwf{n}) / sqrt(length(bwf{n}));
        p75wf(n) = mean(bwf{n}) + std(bwf{n}) / sqrt(length(bwf{n}));
    end
end
udc = udc(1:size(mwf,2));
figure; errorbar(udc, mwf, mwf-p25wf, p75wf-mwf, ':'); hold on
X = [ones(size(udc')) udc'];
[b,bint,r,rint,stats]= regress(  mwf', X);
yfit = b(1)+ b(2)* udc;
plot(udc, yfit, 'r')
clear bwf mwf p25wf p75wf

% Convergence strength distribution
figure
hist(conv_strength_norm(:),24)
[nm0 xout0] = hist(conv_strength_norm(:),24);
[nm xout] = hist(conv_strength_rand_norm(:),xout0);
hold on
plot(xout(nm>0),nm(nm>0),'c')

% Divergence strength distribution
figure
hist(div_strength_norm(:),24)
[nm0 xout0] = hist(div_strength_norm(:),24);
[nm xout] = hist(div_strength_rand_norm(:),xout0);
hold on
plot(xout(nm>0),nm(nm>0),'c')

% Convergence-divergence correlations
figure
plot(conv_strength/60000,div_strength/1200,'k.','MarkerSize',8)
X = [ones(size(conv_strength')) conv_strength'/60000];
[b,bint,r,rint,stats]= regress(div_strength'/1200, X);
yfit = b(1)+ b(2)* conv_strength / 60000;
hold on
plot(conv_strength/60000, yfit, 'r')

figure
plot(conv_strength_norm,div_strength_norm,'k.','MarkerSize',8)
X = [ones(size(conv_strength_norm')) conv_strength_norm'];
[b,bint,r,rint,stats]= regress(div_strength_norm', X);
yfit = b(1)+ b(2)* conv_strength_norm;
hold on
plot(conv_strength_norm, yfit, 'r')

1

% figure;plot3(-conv_strength1/60000,-conv_strength2/60000,conv_strength3/60000,'.')
% grid on
% 
% figure;plot3(-div_strength1/60000,-div_strength2/60000,div_strength3/60000,'.')
% grid on


% -------------------------------------------------------------------------
function [sCnv sDiv convdm divdm wf wt] = lload(pat,patno,eg,nm_rows,nm_cols)

% Directories
global DATADIR
global DATAPATH
inpdir = [DATAPATH 'Ulbert\OITI_' patno '_EEG_' eg '\Contrib2\'];
try
    resdir = [DATAPATH 'Ulbert\OITI_' patno '_EEG_' eg '\Convdivgroup'];
    cd(resdir)
catch
    mkdir(resdir)
    cd(resdir)
end

% Load
fn = [inpdir 'convdiv.mat'];    % MI map
load(fn)

% -------------------------------------------------------------------------
function [sCnv_rand sDiv_rand] = lload2(pat,patno,eg,nm_rows,nm_cols)

% Directories
global DATADIR
global DATAPATH
inpdir = [DATAPATH 'Ulbert\OITI_' patno '_EEG_' eg '\Traject_randconvdiv\'];

% Open (conv)
fn = [inpdir 'rnd_conv_' eg '.fig'];
open(fn)
im = findobj(allchild(gca),'Type','Image');
sCnv_rand = get(im,'Cdata');

% Open (div)
fn = [inpdir 'rnd_div_' eg '.fig'];
open(fn)
im = findobj(allchild(gca),'Type','Image');
sDiv_rand = get(im,'Cdata');