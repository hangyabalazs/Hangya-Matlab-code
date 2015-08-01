function lfburststat
%LFBURSTSTAT   Burst statistics.
%   LFBURSTSTAT calculates summary statistics for burst parameters. Output:
%   burstiness, intraburst frequency, intraburst spike number, burst
%   frequency and firing rate. Box plots are saved with the results of
%   Mann-Whitney tests.
%
%   See also AALL.

% Directories
global DATAPATH
global DATADIR
inpdir1 = [DATAPATH 'LentiFenti\Cluster\mat2_Contr\'];
inpdir2 = [DATAPATH 'LentiFenti\Cluster\mat2_Cre\'];
inpdir_disc1 = [DATADIR 'LentiFenti\disc_Contr\'];
inpdir_disc2 = [DATADIR 'LentiFenti\disc_Cre\'];
resdir = [DATAPATH 'LentiFenti\burststat\'];
mm = pwd;

% Burst statistics
[brstness_Contr ibfr_Contr ibspno_Contr bl_Contr bf_Contr fr_Contr] = main(inpdir1,inpdir_disc1);
[brstness_Cre ibfr_Cre ibspno_Cre bl_Cre bf_Cre fr_Cre] = main(inpdir2,inpdir_disc2);
    
% Plot & save
loc = [];
cd(resdir)
H = burstbox(brstness_Contr,brstness_Cre,loc);
title('Burstiness')
saveas(H,'Burstiness_boxplot.fig');
saveas(H,'Burstiness_boxplot.jpg');
saveas(H,'Burstiness_boxplot.eps');

H = burstbox(ibfr_Contr,ibfr_Cre,loc);
title('IntraBurstFrequency')
saveas(H,'IntraBurstFrequency_boxplot.fig');
saveas(H,'IntraBurstFrequency_boxplot.jpg');
saveas(H,'IntraBurstFrequency_boxplot.eps');

H = burstbox(ibspno_Contr,ibspno_Cre,loc);
title('IntraBurstSpikeNumber')
saveas(H,'IntraBurstSpikeNumber_boxplot.fig');
saveas(H,'IntraBurstSpikeNumber_boxplot.jpg');
saveas(H,'IntraBurstSpikeNumber_boxplot.eps');

H = burstbox(bl_Contr,bl_Cre,loc);
title('BurstLength')
saveas(H,'BurstLength_boxplot.fig');
saveas(H,'BurstLength_boxplot.jpg');
saveas(H,'BurstLength_boxplot.eps');

H = burstbox(bf_Contr,bf_Cre,loc);
title('BurstFrequency')
saveas(H,'BurstFrequency_boxplot.fig');
saveas(H,'BurstFrequency_boxplot.jpg');
saveas(H,'BurstFrequency_boxplot.eps');

H = burstbox(fr_Contr,fr_Cre,loc);
title('FiringRate')
saveas(H,'FiringRate_boxplot.fig');
saveas(H,'FiringRate_boxplot.jpg');
saveas(H,'FiringRate_boxplot.eps');
close all
cd(mm)



% -------------------------------------------------------------------------
function [brstness ibfr ibspno bl bf fr] = main(inpdir,inpdir_disc)

% Filelist
[files files_short] = filelist(inpdir);
sf = length(files);

% MAIN
sr = 20000;         % sampling rate
brstness = nan(1,sf);
ibfr = nan(1,sf);
ibspno = nan(1,sf);
bl = nan(1,sf);
bf = nan(1,sf);
fr = nan(1,sf);
for o = 1:sf
    fname = files_short{o};
    ff3 = [inpdir fname(1:end-4) '_CLUST2'];
    load(ff3)       % load burst analysis results
    ff2 = [inpdir_disc fname(1:end-4) '_d.mat'];
    load(ff2)           % load discriminated unit

    % Burst statistics
    vburst = vdisc(Burst);
    lvb = vburst;
    burstnum = size(lvb,2);
    intraburstiv = [];
    intraburstnum = zeros(1,burstnum);
    for j = 1:burstnum      % computing intraburstiv
        b = vdisc(vdisc>=lvb(1,j)&vdisc<=lvb(2,j));
        db = diff(b);
        intraburstiv = [intraburstiv db];
        intraburstnum(j) = length(b);   % intraburst spike number
    end
    burstiness = (length(intraburstiv) + burstnum) / length(vdisc);
    burstlength = (lvb(2,:) - lvb(1,:)) / sr;
    if ~isempty(intraburstnum)
        intraburstfreq = mean((intraburstnum-1)./burstlength);
        ibspnop = mean(intraburstnum);
    else
        intraburstfreq = NaN;
        ibspnop = NaN;
        warning('lfburststat:NoBursts','Empty variable.')
    end
    burstfreq = sr * (burstnum - 1)  / (vdisc(Burst(2,end)) - vdisc(Burst(1,1)));
    efflen = (vdisc(end) - vdisc(1)) / sr;
    frate = length(vdisc) / efflen;
    
    brstness(o) = burstiness;
    ibfr(o) = intraburstfreq;
    ibspno(o) = ibspnop;
    bl(o) = mean(burstlength) * 1000;    % in ms
    bf(o) = burstfreq;
    fr(o) = frate;
end



% -------------------------------------------------------------------------
function [files2 files2_short] = filelist(inpdir)

% List of filenames
files = dir(inpdir);
files = files(3:end);
files2 = struct('name',[],'date',[],'bytes',[],'isdir',[],'datenum',[]);
files2_short = {};
for i = 1:length(files)
    if ~files(i).isdir
        files2(end+1) = files(i);
        files2_short{end+1} = [files(i).name(1:end-11) '.mat'];
    end
end
files2 = files2(2:end);



% -------------------------------------------------------------------------
function H = burstbox(m1,m2,loc)

H = figure;
boxplot([m1 m2],[zeros(size(m1)) ones(size(m2))],'labels',[{'Control'} {'Cre'}]);
[Wp_eu,Wh_eu] = ranksum(m1,m2,'alpha',0.05);
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
if ~isempty(loc)
    text(tpos1,tpos2,loc)
end