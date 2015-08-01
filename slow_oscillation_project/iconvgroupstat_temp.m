function iconvgroupstat
%ICONCVGROUPSTAT   Group statistics for convergence and divergence maps.
%   ICONVGROUPSTAT calculates convergence and divergence strength
%   distributions as well as distance-contribution strength functions.

% Read Excel file
global DATAPATH
fn = [DATAPATH 'Ulbert\patients.xls'];
headerrows = 0;
[mtx ntx atx] = xlsread(fn);
ntx(1:headerrows,:) = [];
atx(1:headerrows,:) = [];

% Call
conv_strength = [];
div_strength = [];
conv_strength_rand = [];
div_strength_rand = [];
dist_conv = [];
dist_div = [];
where_from = [];
where_to = [];
T = [1:3 7:9 13:15];
T = [4:6 10:12];
T = 1:15;
for k = T
    k
    pat = atx{k,1};
    patno = num2str(atx{k,2});
    eg = num2str(atx{k,3});
    nm_rows = atx{k,4};
    nm_cols = atx{k,5};
    
    [sCnv sDiv convdm divdm wf wt] = lload(pat,patno,eg,nm_rows,nm_cols);
    [sCnv_rand sDiv_rand] = lload2(pat,patno,eg,nm_rows,nm_cols);
    conv_strength = [conv_strength sCnv];
    div_strength = [div_strength sDiv];
    conv_strength_rand = [conv_strength_rand sCnv_rand(1:16)];
    div_strength_rand = [div_strength_rand sDiv_rand(1:16)];
    dist_conv = [dist_conv convdm];
    dist_div = [dist_div divdm];
    where_from = [where_from wf];
    where_to = [where_to wt];
    close all
end

conv_strength = conv_strength(conv_strength>0);
div_strength = div_strength(div_strength>0);
conv_strength_rand = conv_strength_rand(conv_strength_rand>0);
div_strength_rand = div_strength_rand(div_strength_rand>0);

figure
hist(conv_strength(:)/60000,24)
[nm0 xout0] = hist(conv_strength(:)/60000,24);
[nm xout] = hist(conv_strength_rand(:)/60000,xout0);
hold on
plot(xout(nm>0),nm(nm>0),'c')

figure
hist(div_strength(:)/600,25)
[nm0 xout0] = hist(div_strength(:)/600,25);
[nm xout] = hist(div_strength_rand(:)/600,xout0);
hold on
plot(xout(nm>0),nm(nm>0),'c')

figure
plot(dist_conv,where_from,'.')
figure
plot(dist_div,where_to,'.')

% -------------------------------------------------------------------------
function [sCnv sDiv convdm divdm wf wt] = lload(pat,patno,eg,nm_rows,nm_cols)

% Directories
global DATADIR
global DATAPATH
inpdir = [DATAPATH 'Ulbert\OITI_' patno '_EEG_' eg '\Contrib\'];
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