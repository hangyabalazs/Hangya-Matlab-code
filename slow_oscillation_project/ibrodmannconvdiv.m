function ibrodmannconvdiv
%IBRODMANNCONVDIV   Convergence and divergence maps' relation to Brodmann areas.
%   IBRODMANNCONVDIV calculates convergence and divergence strength
%   distributions across Brodmann areas.
%
%   See also ICONVGROUPSTAT2.

% Read Excel file
global DATAPATH
fn = [DATAPATH 'Ulbert\patients_all.xls'];
headerrows = 0;
[mtx ntx atx] = xlsread(fn);
ntx(1:headerrows,:) = [];
atx(1:headerrows,:) = [];
inpdir_Br = [DATAPATH 'Ulbert\reko\'];

% Call
conv_strength = [];
div_strength = [];
conv_strength_norm = [];
div_strength_norm = [];
Brodmann = [];
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
    conv_strength = [conv_strength sCnv(:)'];
    div_strength = [div_strength sDiv(:)'];
    conv_strength_norm = [conv_strength_norm (sCnv(:)'/60000-mean(sCnv_rand(:)/60000))/std(sCnv_rand(:)/60000)];
    div_strength_norm = [div_strength_norm (sDiv(:)'/1200-mean(sDiv_rand(:)/1200))/std(sDiv_rand(:)/1200)];
    close all
    
    fn = [inpdir_Br 'OITI_' patno '_Brodmann.mat'];
    load(fn)
    Br = Br';
    Brodmann = [Brodmann Br(:)'];
end

% I = conv_strength > 0 | div_strength > 0;
% conv_strength = conv_strength(I);
% div_strength = div_strength(I);
% conv_strength_norm = conv_strength_norm(I);
% div_strength_norm = div_strength_norm(I);
% Brodmann = Brodmann(I);

brs = sort(unique(Brodmann));
lenb = length(brs);
conv_Br = zeros(1,lenb);
conv_Br_sd = zeros(1,lenb);
conv_Br_se = zeros(1,lenb);
div_Br = zeros(1,lenb);
div_Br_sd = zeros(1,lenb);
div_Br_se = zeros(1,lenb);
for k = 1:lenb
    inx = Brodmann == brs(k);
    conv_Br(k) = mean(conv_strength_norm(inx));
    conv_Br_sd(k) = std(conv_strength_norm(inx));
    conv_Br_se(k) = std(conv_strength_norm(inx)) / sqrt(length(conv_strength_norm(inx)));
    div_Br(k) = mean(div_strength_norm(inx));
    div_Br_sd(k) = std(div_strength_norm(inx));
    div_Br_se(k) = std(div_strength_norm(inx)) / sqrt(length(div_strength_norm(inx)));
end

figure;plot(conv_Br,'ko')
hold on
errorbar(conv_Br,conv_Br_se,'ko','LineWidth',2)
line([0 25],[0 0],'Color','r','LineWidth',2)

figure;plot(div_Br,'ko')
hold on
errorbar(div_Br,div_Br_se,'ko','LineWidth',2)
line([0 28],[0 0],'Color','r','Linewidth',2)

inxs = [1 2 3 4 5 16 14 15 17 6 7 8 9 10 11 13 12] + 1;
xs = [1 2 3 6 7 8 11 12 13 16 17 20 21 22 23 24 25];
figure;plot(xs,conv_Br(inxs),'ko')
hold on
errorbar(xs,conv_Br(inxs),conv_Br_se(inxs),'ko','LineWidth',2)
line([0 28],[0 0],'Color','r','LineWidth',2)
iinx = [1 2 3 11:17];
errorbar(xs(iinx),conv_Br(inxs(iinx)),conv_Br_se(inxs(iinx)),'o','LineWidth',2,'Color',[0.7 0.7 0.7])
set(gca,'LineWidth',2,'XTick',[])
box off

inxs = [1 2 3 4 5 16 14 15 17 6 7 8 9 10 11 13 12] + 1;
xs = [1 2 3 6 7 8 11 12 13 16 17 20 21 22 23 24 25];
figure;plot(xs,div_Br(inxs),'ko')
hold on
errorbar(xs,div_Br(inxs),div_Br_se(inxs),'ko','LineWidth',2)
line([0 28],[0 0],'Color','r','LineWidth',2)
iinx = [1 2 3 11:17];
errorbar(xs(iinx),div_Br(inxs(iinx)),div_Br_se(inxs(iinx)),'o','LineWidth',2,'Color',[0.7 0.7 0.7])
set(gca,'LineWidth',2,'XTick',[])
box off


% -------------------------------------------------------------------------
function [sCnv sDiv convdm divdm wf wt] = lload(pat,patno,eg,nm_rows,nm_cols)

% Directories
global DATADIR
global DATAPATH
inpdir = [DATAPATH 'Ulbert\OITI_' patno '_EEG_' eg '\Contrib2\'];

% Load
fn = [inpdir 'convdiv.mat'];    % MI map
load(fn)

% -------------------------------------------------------------------------
function [sCnv_rand sDiv_rand] = lload2(pat,patno,eg,nm_rows,nm_cols)

% Directories
global DATADIR
global DATAPATH
inpdir = [DATAPATH 'Ulbert\OITI_' patno '_EEG_' eg '\Traject\'];

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