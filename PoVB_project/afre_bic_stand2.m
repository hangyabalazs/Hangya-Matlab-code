function afre_bic_stand2(inpdir0,atx)
%AFRE_BIC_STAND2   Comparison of baseline and bicuculline.
%   AFRE_BIC_STAND2(DR,ATX) calculates and saves burst parameters and phase
%   statistics for baseline and bicuculline data. Input directory and a
%   matrix containing localization information should be given as an
%   argument (DR, ATX).
%
%   See also APHASE and APHASERUN_BURST2.

% Input argument check
error(nargchk(2,2,nargin))

% Directories
global DATAPATH
inpdir_bas = [inpdir0 'bas\'];
inpdir_bic = [inpdir0 'bic\'];
inpdir1 = [DATAPATH 'Andi\Ketxyl\FreBandRestrict_phase_stand\'];   % phase analysis data
inpdir2 = [DATAPATH 'Andi\Ketxyl\FreBandRestrict_burst_bic_stand\'];   % burst analysis data
resdir = [DATAPATH 'Andi\Ketxyl\AfreBic_stand\'];
mm = pwd;
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
saveas(H,[fname '_Phase.fig']);
saveas(H,[fname '_Phase.eps']);

% Compare burst parameters of bicuculline and baseline
H = burstfig(brstness_bas,brstness_bic,loc);
title('Burstiness')
saveas(H,[fname '_Burstiness.fig']);
saveas(H,[fname '_Burstiness.eps']);

H = burstfig(ibfr_bas,ibfr_bic,loc);
title('IntraBurstFrequency')
saveas(H,[fname '_IntraBurstFrequency.fig']);
saveas(H,[fname '_IntraBurstFrequency.eps']);
H = figure;
subplot(2,1,1)
hist(ibfr_bas,100:20:700)
xlim([100 700])
title('IntraBurstFrequency')
subplot(2,1,2)
hist(ibfr_bic,100:20:700)
xlim([100 700])
title(loc)
saveas(H,[fname '_IntraBurstFrequency_hist.fig']);
saveas(H,[fname '_IntraBurstFrequency_hist.eps']);

H = burstfig(ibspno_bas+rand(size(ibspno_bas))*0.4-0.2,...
    ibspno_bic+rand(size(ibspno_bic))*0.4-0.2,loc);
title('IntraBurstSpikeNumber')
saveas(H,[fname '_IntraBurstSpikeNumber.fig']);
saveas(H,[fname '_IntraBurstSpikeNumber.eps']);
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
saveas(H,[fname '_IntraBurstSpikeNumber_hist.fig']);
saveas(H,[fname '_IntraBurstSpikeNumber_hist.eps']);

H = burstfig(bl_bas,bl_bic,loc);
title('BurstLength')
saveas(H,[fname '_BurstLength.fig']);
saveas(H,[fname '_BurstLength.eps']);
H = figure;
subplot(2,1,1)
hist(bl_bas,0:0.001:0.03)
xlim([0 0.03])
title('BurstLength')
subplot(2,1,2)
hist(bl_bic,0:0.001:0.03)
xlim([0 0.03])
title(loc)
saveas(H,[fname '_BurstLength_hist.fig']);
saveas(H,[fname '_BurstLength_hist.eps']);

H = burstfig(bf_bas,bf_bic,loc);
title('BurstFrequency')
saveas(H,[fname '_BurstFrequency.fig']);
saveas(H,[fname '_BurstFrequency.eps']);

H = burstfig(fr_bas,fr_bic,loc);
title('FiringRate')
saveas(H,[fname '_FiringRate.fig']);
saveas(H,[fname '_FiringRate.eps']);

H = burstbox(brstness_bas,brstness_bic,loc);
title('Burstiness')
saveas(H,[fname '_Burstiness_boxplot.fig']);
saveas(H,[fname '_Burstiness_boxplot.eps']);

H = burstbox(ibfr_bas,ibfr_bic,loc);
title('IntraBurstFrequency')
saveas(H,[fname '_IntraBurstFrequency_boxplot.fig']);
saveas(H,[fname '_IntraBurstFrequency_boxplot.eps']);

H = burstbox(bl_bas,bl_bic,loc);
title('BurstLength')
saveas(H,[fname '_BurstLength_boxplot.fig']);
saveas(H,[fname '_BurstLength_boxplot.eps']);

H = burstbox(bf_bas,bf_bic,loc);
title('BurstFrequency')
saveas(H,[fname '_BurstFrequency_boxplot.fig']);
saveas(H,[fname '_BurstFrequency_boxplot.eps']);

H = burstbox(fr_bas,fr_bic,loc);
title('FiringRate')
saveas(H,[fname '_FiringRate_boxplot.fig']);
saveas(H,[fname '_FiringRate_boxplot.eps']);
close all
cd(mm)

% -------------------------------------------------------------------------
function [aang_afsp, meanrl, Burstiness, ibfr, ibspno, bl, BurstFrequency,...
    FiringRate fname] = impdata(inpdir1,inpdir2,files_short,sf)

% Load
fname = files_short{sf}     % filename
fname = fname(1:end-4);
ff1 = [inpdir1 fname '_PHASE.mat'];
load(ff1)
ff2 = [inpdir2 fname '_CLUST2.mat'];
load(ff2)

% Mean resultant length
[ftm mn meanrl] = mvlmn(aang_afsp,'deg');

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