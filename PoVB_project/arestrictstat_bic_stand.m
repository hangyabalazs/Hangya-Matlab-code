function arestrictstat_bic_stand
%ARESTRICTSTAT_BIC_STAND   Comparison of baseline and bicuculline.
%   ARESTRICTSTAT_BIC_STAND calculates and saves burst parameters and phase
%   statistics for baseline and bicuculline data. Cells are grouped based
%   on their anatomical location.
%
%   See also AFRE_RESTRICT_BIC_STAND2.

% Directories
tabledir = ['Y:\_Projects\AUJ_ISTVAN\TABLES\'];
inproot = 'Y:\_Projects\AUJ_ISTVAN\DATA\MAT\mat_ket_xyl\';
dr = dir(inproot);
inpadd = {};
for k = 3:length(dr)
    if dr(k).isdir
        inpadd{end+1} = dr(k).name;
    end
end

% Read Excel file
fn = [tabledir 'tablazat_balazsnak.xls'];
headerrows = 0;
[mtx ntx atx] = xlsread(fn);
ntx(1:headerrows,:) = [];
atx(1:headerrows,:) = [];

% Call main
AANG_BAS = struct('Po',{[]},'VB',{[]},'PoVPM',{[]},'LD',{[]},'nRT',{[]});
MEANRL_BAS = struct('Po',{[]},'VB',{[]},'PoVPM',{[]},'LD',{[]},'nRT',{[]});
BRSTNESS_BAS = struct('Po',{[]},'VB',{[]},'PoVPM',{[]},'LD',{[]},'nRT',{[]});
IBFR_BAS = struct('Po',{[]},'VB',{[]},'PoVPM',{[]},'LD',{[]},'nRT',{[]});
IBSPNO_BAS = struct('Po',{[]},'VB',{[]},'PoVPM',{[]},'LD',{[]},'nRT',{[]});
BL_BAS = struct('Po',{[]},'VB',{[]},'PoVPM',{[]},'LD',{[]},'nRT',{[]});
BF_BAS = struct('Po',{[]},'VB',{[]},'PoVPM',{[]},'LD',{[]},'nRT',{[]});
FR_BAS = struct('Po',{[]},'VB',{[]},'PoVPM',{[]},'LD',{[]},'nRT',{[]});
AANG_BIC = struct('Po',{[]},'VB',{[]},'PoVPM',{[]},'LD',{[]},'nRT',{[]});
MEANRL_BIC = struct('Po',{[]},'VB',{[]},'PoVPM',{[]},'LD',{[]},'nRT',{[]});
BRSTNESS_BIC = struct('Po',{[]},'VB',{[]},'PoVPM',{[]},'LD',{[]},'nRT',{[]});
IBFR_BIC = struct('Po',{[]},'VB',{[]},'PoVPM',{[]},'LD',{[]},'nRT',{[]});
IBSPNO_BIC = struct('Po',{[]},'VB',{[]},'PoVPM',{[]},'LD',{[]},'nRT',{[]});
BL_BIC = struct('Po',{[]},'VB',{[]},'PoVPM',{[]},'LD',{[]},'nRT',{[]});
BF_BIC = struct('Po',{[]},'VB',{[]},'PoVPM',{[]},'LD',{[]},'nRT',{[]});
FR_BIC = struct('Po',{[]},'VB',{[]},'PoVPM',{[]},'LD',{[]},'nRT',{[]});
for k = 1:length(inpadd)
    inpdir = [inproot inpadd{k} '\']
    [AANG_BAS, MEANRL_BAS, BRSTNESS_BAS, IBFR_BAS, IBSPNO_BAS, BL_BAS, BF_BAS, FR_BAS...
        AANG_BIC, MEANRL_BIC, BRSTNESS_BIC, IBFR_BIC, IBSPNO_BIC, BL_BIC, BF_BIC, FR_BIC] ...
        = main(inpdir,atx,AANG_BAS,MEANRL_BAS,BRSTNESS_BAS,IBFR_BAS,IBSPNO_BAS,BL_BAS,BF_BAS,FR_BAS,...
        AANG_BIC,MEANRL_BIC,BRSTNESS_BIC,IBFR_BIC,IBSPNO_BIC,BL_BIC,BF_BIC,FR_BIC)
   
    
end
basbic_compare(AANG_BAS.Po, MEANRL_BAS.Po, BRSTNESS_BAS.Po, IBFR_BAS.Po, IBSPNO_BAS.Po, BL_BAS.Po, BF_BAS.Po, FR_BAS.Po,...
    AANG_BIC.Po, MEANRL_BIC.Po, BRSTNESS_BIC.Po, IBFR_BIC.Po, IBSPNO_BIC.Po, BL_BIC.Po, BF_BIC.Po, FR_BIC.Po, 'Po')
basbic_compare(AANG_BAS.VB, MEANRL_BAS.VB, BRSTNESS_BAS.VB, IBFR_BAS.VB, IBSPNO_BAS.VB, BL_BAS.VB, BF_BAS.VB, FR_BAS.VB,...
    AANG_BIC.VB, MEANRL_BIC.VB, BRSTNESS_BIC.VB, IBFR_BIC.VB, IBSPNO_BIC.VB, BL_BIC.VB, BF_BIC.VB, FR_BIC.VB, 'VB')
basbic_compare(AANG_BAS.PoVPM, MEANRL_BAS.PoVPM, BRSTNESS_BAS.PoVPM, IBFR_BAS.PoVPM, IBSPNO_BAS.PoVPM, BL_BAS.PoVPM, BF_BAS.PoVPM, FR_BAS.PoVPM,...
    AANG_BIC.PoVPM, MEANRL_BIC.PoVPM, BRSTNESS_BIC.PoVPM, IBFR_BIC.PoVPM, IBSPNO_BIC.PoVPM, BL_BIC.PoVPM, BF_BIC.PoVPM, FR_BIC.PoVPM, 'PoVPM')
basbic_compare(AANG_BAS.LD, MEANRL_BAS.LD, BRSTNESS_BAS.LD, IBFR_BAS.LD, IBSPNO_BAS.LD, BL_BAS.LD, BF_BAS.LD, FR_BAS.LD,...
    AANG_BIC.LD, MEANRL_BIC.LD, BRSTNESS_BIC.LD, IBFR_BIC.LD, IBSPNO_BIC.LD, BL_BIC.LD, BF_BIC.LD, FR_BIC.LD, 'LD')
basbic_compare(AANG_BAS.nRT, MEANRL_BAS.nRT, BRSTNESS_BAS.nRT, IBFR_BAS.nRT, IBSPNO_BAS.nRT, BL_BAS.nRT, BF_BAS.nRT, FR_BAS.nRT,...
    AANG_BIC.nRT, MEANRL_BIC.nRT, BRSTNESS_BIC.nRT, IBFR_BIC.nRT, IBSPNO_BIC.nRT, BL_BIC.nRT, BF_BIC.nRT, FR_BIC.nRT, 'nRT')

% -------------------------------------------------------------------------
function [AANG_BAS, MEANRL_BAS, BRSTNESS_BAS, IBFR_BAS, IBSPNO_BAS, BL_BAS, BF_BAS, FR_BAS...
    AANG_BIC, MEANRL_BIC, BRSTNESS_BIC, IBFR_BIC, IBSPNO_BIC, BL_BIC, BF_BIC, FR_BIC] ...
    = main(inpdir0,atx,AANG_BAS,MEANRL_BAS,BRSTNESS_BAS,IBFR_BAS,IBSPNO_BAS,BL_BAS,BF_BAS,FR_BAS,...
    AANG_BIC,MEANRL_BIC,BRSTNESS_BIC,IBFR_BIC,IBSPNO_BIC,BL_BIC,BF_BIC,FR_BIC)

% Directories
global DATAPATH
inpdir_bas = [inpdir0 'bas\'];
inpdir_bic = [inpdir0 'bic\'];
inpdir1 = [DATAPATH 'Andi\Ketxyl\FreBandRestrict_phase_stand\'];   % phase analysis data
inpdir2 = [DATAPATH 'Andi\Ketxyl\FreBandRestrict_burst_bic_stand\'];   % burst analysis data
dbstop if error

% Filelist
[files_bas files_short_bas] = filelist(inpdir_bas);
[files_bic files_short_bic] = filelist(inpdir_bic);
sf_bas = length(files_short_bas);
sf_bic = length(files_short_bic);
if isequal(sf_bas,0)
    disp('No bas file.')
    return
end
if isequal(sf_bic,0)
    disp('No bic file.')
    return
end
fname = files_short_bas{sf_bas};
fname = fname(1:end-4);
ff1 = [inpdir1 fname '_PHASE.mat'];
if ~exist(ff1,'file')
    disp('No phase data.')
    return
end
ff2 = [inpdir2 fname '_CLUST2.mat'];
if ~exist(ff2,'file')
    disp('No burst data.')
    return
end
fname = files_short_bic{sf_bic};
fname = fname(1:end-4);
ff1 = [inpdir1 fname '_PHASE.mat'];
if ~exist(ff1,'file')
    disp('No phase data.')
    return
end
ff2 = [inpdir2 fname '_CLUST2.mat'];
if ~exist(ff2,'file')
    disp('No burst data.')
    return
end

% Import
[aang_bas, meanrl_bas, brstness_bas, ibfr_bas, ibspno_bas, bl_bas, bf_bas, fr_bas]...
    = impdata(inpdir1,inpdir2,files_short_bas,sf_bas);
[aang_bic, meanrl_bic, brstness_bic, ibfr_bic, ibspno_bic, bl_bic, bf_bic, fr_bic fname]...
    = impdata(inpdir1,inpdir2,files_short_bic,sf_bic);

% Localization
cmps = strread(fname,'%s','delimiter','_.');
fname2 = [cmps{1} '_' cmps{2}];
inx = [find(strcmp({atx{:,1}},fname2)) find(strcmp({atx{:,2}},fname2))];
loc = atx{inx,3};

% Create pooled samples
dsc1 = 100;
dsc2 = 10;
if length(aang_bas) >= dsc1 && length(aang_bic) >= dsc1
    rp = randperm(length(aang_bas));     % draw random sample for balanced pooled sample
    aang_bas_ds = aang_bas(rp(1:dsc1));
    eval(['AANG_BAS.' loc ' = [AANG_BAS.' loc ' aang_bas_ds];']);
    rp = randperm(length(aang_bic));
    aang_bic_ds = aang_bic(rp(1:dsc1));
    eval(['AANG_BIC.' loc ' = [AANG_BIC.' loc ' aang_bic_ds];']);
end
if length(meanrl_bas) >= dsc2 && length(meanrl_bic) >= dsc2
    rp = randperm(length(meanrl_bas));
    meanrl_bas_ds = meanrl_bas(rp(1:dsc2));
    eval(['MEANRL_BAS.' loc ' = [MEANRL_BAS.' loc ' meanrl_bas_ds];']);
    rp = randperm(length(meanrl_bic));
    meanrl_bic_ds = meanrl_bic(rp(1:dsc2));
    eval(['MEANRL_BIC.' loc ' = [MEANRL_BIC.' loc ' meanrl_bic_ds];']);
end
if length(brstness_bas) >= dsc2 && length(brstness_bic) >= dsc2
    rp = randperm(length(brstness_bas));
    brstness_bas_ds = brstness_bas(rp(1:dsc2));
    eval(['BRSTNESS_BAS.' loc ' = [BRSTNESS_BAS.' loc ' brstness_bas_ds];']);
    rp = randperm(length(brstness_bic));
    brstness_bic_ds = brstness_bic(rp(1:dsc2));
    eval(['BRSTNESS_BIC.' loc ' = [BRSTNESS_BIC.' loc ' brstness_bic_ds];']);
end
if length(ibfr_bas) >= dsc1 && length(ibfr_bic) >= dsc1
    rp = randperm(length(ibfr_bas));
    ibfr_bas_ds = ibfr_bas(rp(1:dsc1));
    eval(['IBFR_BAS.' loc ' = [IBFR_BAS.' loc ' ibfr_bas_ds];']);
    rp = randperm(length(ibfr_bic));
    ibfr_bic_ds = ibfr_bic(rp(1:dsc1));
    eval(['IBFR_BIC.' loc ' = [IBFR_BIC.' loc ' ibfr_bic_ds];']);
end
if length(ibspno_bas) >= dsc1 && length(ibspno_bic) >= dsc1
    rp = randperm(length(ibspno_bas));
    ibspno_bas_ds = ibspno_bas(rp(1:dsc1));
    eval(['IBSPNO_BAS.' loc ' = [IBSPNO_BAS.' loc ' ibspno_bas_ds];']);
    rp = randperm(length(ibspno_bic));
    ibspno_bic_ds = ibspno_bic(rp(1:dsc1));
    eval(['IBSPNO_BIC.' loc ' = [IBSPNO_BIC.' loc ' ibspno_bic_ds];']);
end
if length(bl_bas) >= dsc1 && length(bl_bic) >= dsc1
    rp = randperm(length(bl_bas));
    bl_bas_ds = bl_bas(rp(1:dsc1));
    eval(['BL_BAS.' loc ' = [BL_BAS.' loc ' bl_bas_ds];']);
    rp = randperm(length(bl_bic));
    bl_bic_ds = bl_bic(rp(1:dsc1));
    eval(['BL_BIC.' loc ' = [BL_BIC.' loc ' bl_bic_ds];']);
end
if length(bf_bas) >= dsc2 && length(bf_bic) >= dsc2
    rp = randperm(length(bf_bas));
    bf_bas_ds = bf_bas(rp(1:dsc2));
    eval(['BF_BAS.' loc ' = [BF_BAS.' loc ' bf_bas_ds];']);
    rp = randperm(length(bf_bic));
    bf_bic_ds = bf_bic(rp(1:dsc2));
    eval(['BF_BIC.' loc ' = [BF_BIC.' loc ' bf_bic_ds];']);
end
if length(fr_bas) >= dsc2 && length(fr_bic) >= dsc2
    rp = randperm(length(fr_bas));
    fr_bas_ds = fr_bas(rp(1:dsc2));
    eval(['FR_BAS.' loc ' = [FR_BAS.' loc ' fr_bas_ds];']);
    rp = randperm(length(fr_bic));
    fr_bic_ds = fr_bic(rp(1:dsc2));
    eval(['FR_BIC.' loc ' = [FR_BIC.' loc ' fr_bic_ds];']);
end

% -------------------------------------------------------------------------
function basbic_compare(aang_bas, meanrl_bas, brstness_bas, ibfr_bas, ibspno_bas, bl_bas, bf_bas, fr_bas,...
    aang_bic, meanrl_bic, brstness_bic, ibfr_bic, ibspno_bic, bl_bic, bf_bic, fr_bic, loc)

% Directories
global DATAPATH
resdir = [DATAPATH 'Andi\Ketxyl\RestrictStat_bic_stand\'];
mm = pwd;

% Compare phase distributions of bicuculline and baseline
edges = -180:20:180;     % edges for phase histogram
cnts = (edges(1:end-1) + edges(2:end)) / 2;
[nm_bas,xout_bas] = histc(aang_bas,edges);   % phase histogram
nm_bas = nm_bas(1:end-1);
[nm_bic,xout_bic] = histc(aang_bic,edges);
nm_bic = nm_bic(1:end-1);
H = figure;
subplot(2,2,1)
B = bar(cnts,nm_bas'/length(aang_bas));
set(B,'FaceColor',[0.16 0.38 0.27])
ylim([0 0.3])
title('bas')
subplot(2,2,3)
B = bar(cnts,nm_bic'/length(aang_bic));
set(B,'FaceColor',[0.16 0.38 0.27])
ylim([0 0.3])
title('bic')
[Wp mn1 mn2 mvl1 mvl2] = wtest(aang_bas'*pi/180,aang_bic'*pi/180);
Kp = kappacompare2(aang_bas,aang_bic,'deg');
str = {['Watson test: ' num2str(Wp(1)) ' < p < ' num2str(Wp(2))]...
    ['bas mean: ' num2str(mn1)] ['bic mean: ' num2str(mn2)]...
    ['bas mvl: ' num2str(mvl1)] ['bic mvl: ' num2str(mvl2)]...
    ['Randomization test: p = ' num2str(Kp)] ['Localization: ' loc]};
uicontrol('Style','text','Unit','normalized','Position',...
    [0.52 0.37 0.4 0.55],'FontSize',12,'HorizontalAlignment',...
    'left','String',str,'BackgroundColor',[1 1 1]);
cd(resdir)
dbclear if error
saveas(H,[loc '_Phase.fig']);
saveas(H,[loc '_Phase.eps']);

% Compare burst parameters of bicuculline and baseline
H = burstfig(brstness_bas,brstness_bic,loc);
title('Burstiness')
saveas(H,[loc '_Burstiness.fig']);
saveas(H,[loc '_Burstiness.eps']);

H = burstfig(ibfr_bas,ibfr_bic,loc);
title('IntraBurstFrequency')
saveas(H,[loc '_IntraBurstFrequency.fig']);
saveas(H,[loc '_IntraBurstFrequency.eps']);
H = figure;
subplot(2,1,1)
hist(ibfr_bas,100:20:700)
xlim([100 700])
title('IntraBurstFrequency')
subplot(2,1,2)
hist(ibfr_bic,100:20:700)
xlim([100 700])
title(loc)
saveas(H,[loc '_IntraBurstFrequency_hist.fig']);
saveas(H,[loc '_IntraBurstFrequency_hist.eps']);

H = burstfig(ibspno_bas+rand(size(ibspno_bas))*0.4-0.2,...
    ibspno_bic+rand(size(ibspno_bic))*0.4-0.2,loc);
title('IntraBurstSpikeNumber')
saveas(H,[loc '_IntraBurstSpikeNumber.fig']);
saveas(H,[loc '_IntraBurstSpikeNumber.eps']);
H = figure;
subplot(2,1,1)
hist(ibspno_bas,1:1:10)
xlim([1 11])
title('IntraBurstSpikeNumber')
subplot(2,1,2)
hist(ibspno_bic,1:1:10)
xlim([1 11])
title(loc)
nmibspno_bas = hist(ibspno_bas,2:1:10);
nmibspno_bic = hist(ibspno_bic,2:1:10);
[Ch,Cp] = b_chi2test2(nmibspno_bas,nmibspno_bic);
y_lim = ylim;
x_lim = xlim;
tpos1 = x_lim(1) + (x_lim(2) - x_lim(1)) * 0.5;
tpos2 = y_lim(1) + (y_lim(2)-y_lim(1)) * 0.8;
text(tpos1,tpos2,['Chi-square test: p = ' num2str(Cp)])
tpos1 = x_lim(1) + (x_lim(2) - x_lim(1)) * 0.5;
tpos2 = y_lim(1) + (y_lim(2)-y_lim(1)) * 0.7;
text(tpos1,tpos2,['Means (bas, bic): ' num2str(mean(ibspno_bas)) ...
    ' ' num2str(mean(ibspno_bic))])
saveas(H,[loc '_IntraBurstSpikeNumber_hist.fig']);
saveas(H,[loc '_IntraBurstSpikeNumber_hist.eps']);

H = burstfig(bl_bas,bl_bic,loc);
title('BurstLength')
saveas(H,[loc '_BurstLength.fig']);
saveas(H,[loc '_BurstLength.eps']);
H = figure;
subplot(2,1,1)
hist(bl_bas,0:0.001:0.03)
xlim([0 0.03])
title('BurstLength')
subplot(2,1,2)
hist(bl_bic,0:0.001:0.03)
xlim([0 0.03])
title(loc)
saveas(H,[loc '_BurstLength_hist.fig']);
saveas(H,[loc '_BurstLength_hist.eps']);

H = burstfig(bf_bas,bf_bic,loc);
title('BurstFrequency')
saveas(H,[loc '_BurstFrequency.fig']);
saveas(H,[loc '_BurstFrequency.eps']);

H = burstfig(fr_bas,fr_bic,loc);
title('FiringRate')
saveas(H,[loc '_FiringRate.fig']);
saveas(H,[loc '_FiringRate.eps']);

if ~isempty(brstness_bas) && ~isempty(brstness_bic)
    H = burstbox(brstness_bas,brstness_bic,loc);
    title('Burstiness')
    saveas(H,[loc '_Burstiness_boxplot.fig']);
    saveas(H,[loc '_Burstiness_boxplot.eps']);
end

if ~isempty(ibfr_bas) && ~isempty(ibfr_bic)
    H = burstbox(ibfr_bas,ibfr_bic,loc);
    title('IntraBurstFrequency')
    saveas(H,[loc '_IntraBurstFrequency_boxplot.fig']);
    saveas(H,[loc '_IntraBurstFrequency_boxplot.eps']);
end

if ~isempty(bl_bas) && ~isempty(bl_bic)
    H = burstbox(bl_bas,bl_bic,loc);
    title('BurstLength')
    saveas(H,[loc '_BurstLength_boxplot.fig']);
    saveas(H,[loc '_BurstLength_boxplot.eps']);
end

if ~isempty(bf_bas) && ~isempty(bf_bic)
    H = burstbox(bf_bas,bf_bic,loc);
    title('BurstFrequency')
    saveas(H,[loc '_BurstFrequency_boxplot.fig']);
    saveas(H,[loc '_BurstFrequency_boxplot.eps']);
end

if ~isempty(fr_bas) && ~isempty(fr_bic)
    H = burstbox(fr_bas,fr_bic,loc);
    title('FiringRate')
    saveas(H,[loc '_FiringRate_boxplot.fig']);
    saveas(H,[loc '_FiringRate_boxplot.eps']);
end
close all
cd(mm)

% -------------------------------------------------------------------------
function [aang_afsp, r_afsp, Burstiness, ibfr, ibspno, bl, BurstFrequency,...
    FiringRate fname] = impdata(inpdir1,inpdir2,files_short,sf)

% Load
fname = files_short{sf}     % filename
fname = fname(1:end-4);
ff1 = [inpdir1 fname '_PHASE.mat'];
load(ff1)
ff2 = [inpdir2 fname '_CLUST2.mat'];
load(ff2)

% Mean resultant length
% [ftm mn meanrl] = mvlmn(aang_afsp,'deg');

% Burst parameters
ibfr = IntraBurstFrequency.all;
ibspno = IntraBurstSpikeNumber.all;
bl = BurstLength.all;

% -------------------------------------------------------------------------
function [files2 files2_short] = filelist(inpdir)

% List of filenames
files = dir(inpdir);
files = files(3:end);
files2 = struct([]);
files2_short = {};
for i = 1:length(files)
    if ~files(i).isdir
        if isempty(files2)
            files2 = files(i);
        else
            files2(end+1) = files(i);
        end
        files2_short{end+1} = files(i).name;
    end
end

% -------------------------------------------------------------------------
function H = burstfig(m1,m2,loc)

H = figure;
hold on
plot(ones(size(m1))+rand(size(m1))*0.6-0.3,m1,'.')
plot(2*ones(size(m2))+rand(size(m2))*0.6-0.3,m2,'.')
xlim([0.5 2.5])
set(gca,'XTick',1:2)
set(gca,'XTickLabel',{'bas','bic'})
y_lim = ylim;
x_lim = xlim;
tpos1 = x_lim(1) + (x_lim(2) - x_lim(1)) / 2;
tpos2 = y_lim(1) + (y_lim(2)-y_lim(1)) * 4 / 5;
text(tpos1,tpos2,loc)

% -------------------------------------------------------------------------
function H = burstbox(m1,m2,loc)

H = figure;
boxplot([m1 m2],[zeros(size(m1)) ones(size(m2))],'labels',[{'bas'} {'bic'}]);
[Wp_eu,Wh_eu] = b_ranksum2(m1,m2,'alpha',0.05);
if Wh_eu
    clr = 'red';
else
    clr = 'black';
end
y_lim = ylim;
x_lim = xlim;
tpos1 = x_lim(1) + (x_lim(2) - x_lim(1)) / 2;
tpos2 = y_lim(1) + (y_lim(2)-y_lim(1)) * 3 / 5;
text(tpos1,tpos2,num2str(Wp_eu),'Color',clr,'Horizontalalignment','center')
tpos1 = x_lim(1) + (x_lim(2) - x_lim(1)) / 2;
tpos2 = y_lim(1) + (y_lim(2)-y_lim(1)) * 4 / 5;
text(tpos1,tpos2,loc)

% -------------------------------------------------------------------------
function [p mn1 mn2 mvl1 mvl2] = wtest(fr1_all,fr2_all)

[U2 p] = b_watsontwo(fr1_all,fr2_all);
ftm = sum(exp(1).^(i*fr1_all)) / length(fr1_all);    % first trigonometric moment
mn1 = angle(ftm);   % mean angle
mn1 = mod(mn1,2*pi) * 180 / pi;
mn1(mn1>180) = mn1 - 360;
mvl1 = abs(ftm);     % mean resultant length
ftm = sum(exp(1).^(i*fr2_all)) / length(fr2_all);    % first trigonometric moment
mn2 = angle(ftm);   % mean angle
mn2 = mod(mn2,2*pi) * 180 / pi;
mn2(mn2>180) = mn2 - 360;
mvl2 = abs(ftm);