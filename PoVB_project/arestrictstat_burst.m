function arestrictstat_burst
%ARESTRICTSTAT_BURST   Burst statistics for bounded frequency band.
%   ARESTRICTSTAT_BURST calculates summary burst statistics for Po, VPM,
%   PoVPM, LD and nRT data. It plots and saves result figures: one
%   displaying mean values of burst parameters for all cell and another
%   representing burst parameters as boxplots. Ordinal interspike interval
%   statistics are also calculated and saved Edit code to modify input and
%   output directories!
%
%   See also ARESTRICTSTAT_PHASE and ARESTRICTSTAT_BURSTPHASE.

% Input argument check
error(nargchk(0,0,nargin))
dbstop if error

% Directories
global DATAPATH
inpdir = 'Y:\_Projects\AUJ_ISTVAN\DATA\MAT\mat_ket_xyl\';
inpdir2 = [DATAPATH 'Andi\Ketxyl\FreBandRestrict_burst_stand\'];   % burst analysis data
tabledir = ['Y:\_Projects\AUJ_ISTVAN\TABLES\'];
resdir = [DATAPATH 'Andi\Ketxyl\RestrictStat_stand\'];
mm = pwd;
dr = dir(inpdir);

% Main
BN = struct('Po',[],'VB',[],'PoVPM',[],'LD',[],'nRT',[]);
IBFR = struct('Po',[],'VB',[],'PoVPM',[],'LD',[],'nRT',[]);
IBSPNO = struct('Po',[],'VB',[],'PoVPM',[],'LD',[],'nRT',[]);
BL = struct('Po',[],'VB',[],'PoVPM',[],'LD',[],'nRT',[]);
BF = struct('Po',[],'VB',[],'PoVPM',[],'LD',[],'nRT',[]);
FR = struct('Po',[],'VB',[],'PoVPM',[],'LD',[],'nRT',[]);
OISI = struct('Po',{{}},'VB',{{}},'PoVPM',{{}},'LD',{{}},'nRT',{{}});
for o = 3:length(dr)
    inpadd = dr(o).name;    % load burst data
    cd(inpdir)
    cd(inpadd)
    cd('bas')
    ddr = dir(pwd);
    fn = ddr(end).name(1:end-4);
    cmps = strread(fn,'%s','delimiter','_');
    fname = [cmps{1} '_' cmps{2}];
    ff = [inpdir2 fn '_CLUST2.mat'];
    try
        load(ff)
    catch
        lasterr
        continue
    end
    ff = [tabledir 'tablazat_Balazsnak'];   % load position data
    [tbl0 tbl] = xlsread(ff);
    inx = find(strcmp({tbl{:,1}},fname));
    loc = tbl{inx,3};
    
    eval(['BN.' loc '(end+1) = Burstiness;']);  % burst parameters
    eval(['IBFR.' loc '(end+1) = IntraBurstFrequency.mean;']);
    eval(['IBSPNO.' loc '(end+1) = IntraBurstSpikeNumber.mean;']);
    eval(['BL.' loc '(end+1) = BurstLength.mean;']);
    eval(['BF.' loc '(end+1) = BurstFrequency;']);
    eval(['FR.' loc '(end+1) = FiringRate;']);
    for k = 1:length(IsiMatrix.mean)
        nt = eval(['length(OISI.' loc ');']);
        if k > nt
            ntt = 0;
        else
            ntt = eval(['size(OISI.' loc '{k},2);']);
        end
        for t = 1:k
            eval(['OISI.' loc '{k}(t,ntt+1) = IsiMatrix.mean(k,t);']);
        end
    end
end

% Plot and save
dbclear if error
cd(resdir)
H = OISIfig(OISI,'Po');
saveas(H,'OrdinalIsi_Po.fig')
saveas(H,'OrdinalIsi_Po.eps')
H = OISIfig(OISI,'VB');
saveas(H,'OrdinalIsi_VB.fig')
saveas(H,'OrdinalIsi_VB.eps')
H = OISIfig(OISI,'PoVPM');
saveas(H,'OrdinalIsi_PoVPM.fig')
saveas(H,'OrdinalIsi_PoVPM.eps')
H = OISIfig(OISI,'LD');
saveas(H,'OrdinalIsi_LD.fig')
saveas(H,'OrdinalIsi_LD.eps')
H = OISIfig(OISI,'nRT');
saveas(H,'OrdinalIsi_nRT.fig')
saveas(H,'OrdinalIsi_nRT.eps')

H = burstfig(BN);
title('Burstiness')
saveas(H,'Burstiness.fig');
saveas(H,'Burstiness.eps');

H = burstfig(IBFR);
title('IntraBurstFrequency')
saveas(H,'IntraBurstFrequency.fig');
saveas(H,'IntraBurstFrequency.eps');

H = burstfig(IBSPNO);
title('IntraBurstSpikeNumber')
saveas(H,'IntraBurstSpikeNumber.fig');
saveas(H,'IntraBurstSpikeNumber.eps');

H = burstfig(BL);
title('BurstLength')
saveas(H,'BurstLength.fig');
saveas(H,'BurstLength.eps');

H = burstfig(BF);
title('BurstFrequency')
saveas(H,'BurstFrequency.fig');
saveas(H,'BurstFrequency.eps');

H = burstfig(FR);
title('FiringRate')
saveas(H,'FiringRate.fig');
saveas(H,'FiringRate.eps');

H = burstbox(BN);
title('Burstiness')
saveas(H,'Burstiness_boxplot.fig');
saveas(H,'Burstiness_boxplot.eps');

H = burstbox(IBFR);
title('IntraBurstFrequency')
saveas(H,'IntraBurstFrequency_boxplot.fig');
saveas(H,'IntraBurstFrequency_boxplot.eps');

H = burstbox(IBSPNO);
title('IntraBurstSpikeNumber')
saveas(H,'IntraBurstSpikeNumber_boxplot.fig');
saveas(H,'IntraBurstSpikeNumber_boxplot.eps');

H = burstbox(BL);
title('BurstLength')
saveas(H,'BurstLength_boxplot.fig');
saveas(H,'BurstLength_boxplot.eps');

H = burstbox(BF);
title('BurstFrequency')
saveas(H,'BurstFrequency_boxplot.fig');
saveas(H,'BurstFrequency_boxplot.eps');

H = burstbox(FR);
title('FiringRate')
saveas(H,'FiringRate_boxplot.fig');
saveas(H,'FiringRate_boxplot.eps');

burststatistics.BN = burststat(BN);
burststatistics.IBFR = burststat(IBFR);
burststatistics.IBSPNO = burststat(IBSPNO);
burststatistics.BL = burststat(BL);
burststatistics.BF = burststat(BF);
burststatistics.FR = burststat(FR);
save burststatistics burststatistics

cd(mm)

% -------------------------------------------------------------------------
function H = burstfig(B)

m1 = B.Po;
m2 = B.VB;
m3 = B.PoVPM;
m4 = B.LD;
m5 = B.nRT;

H = figure;
hold on
plot(ones(size(m1)),m1,'.')
plot(2*ones(size(m2)),m2,'.')
plot(3*ones(size(m3)),m3,'.')
plot(4*ones(size(m4)),m4,'.')
plot(5*ones(size(m5)),m5,'.')
xlim([0 6])
set(gca,'XTick',1:5)
set(gca,'XTickLabel',{'Po','VB','PoVPM','LD','nRT'})

% -------------------------------------------------------------------------
function H = burstbox(B)

m1 = B.Po;
m2 = B.VB;
m3 = B.PoVPM;
m4 = B.LD;
m5 = B.nRT;

H = figure;
boxplot([m1 m2 m3 m4 m5],[zeros(size(m1)) ones(size(m2)) 2*ones(size(m3))...
    3*ones(size(m4)) 4*ones(size(m5))],'labels',[{'Po'} {'VB'} {'PoVPM'}...
    {'LD'} {'nRT'}]);
[Wp_eu,Wh_eu] = b_ranksum2(m1,m2,'alpha',0.05);
if Wh_eu
    clr = 'red';
else
    clr = 'black';
end
y_lim = ylim;
x_lim = xlim;
tpos1 = x_lim(1) + (x_lim(2) - x_lim(1)) / 5;
tpos2 = y_lim(1) + (y_lim(2)-y_lim(1)) * 4 / 5;
text(tpos1,tpos2,num2str(Wp_eu),'Color',clr,'Horizontalalignment','center')

% -------------------------------------------------------------------------
function S = burststat(B)

m1 = B.Po;
m2 = B.VB;
m3 = B.PoVPM;
m4 = B.LD;
m5 = B.nRT;

S.median_Po = median(m1);
S.median_VB = median(m2);
S.median_PoVPM = median(m3);
S.median_LD = median(m4);
S.median_nRT = median(m5);
[S.p_Po_VB] = b_ranksum2(m1,m2,'alpha',0.05);
[S.p_Po_PoVPM] = b_ranksum2(m1,m3,'alpha',0.05);
[S.p_Po_LD] = b_ranksum2(m1,m4,'alpha',0.05);
[S.p_VB_PoVPM] = b_ranksum2(m2,m3,'alpha',0.05);
[S.p_VB_LD] = b_ranksum2(m2,m4,'alpha',0.05);
[S.p_PoVPM_LD] = b_ranksum2(m3,m4,'alpha',0.05);

% -------------------------------------------------------------------------
function H = OISIfig(OISI,Bstr)

H = figure;
B = eval(['OISI.' Bstr]);
hold on
for t = 1:length(B);
    mk = size(B{t},1);
    mn = zeros(1,mk);
    SE = zeros(1,mk);
    for k = 1:mk
        mn(k) = mean(B{t}(k,:));
        SE(k) = std(B{t}(k,:)) / sqrt(size(B{t},2));
    end
    plot(mn,'k')
    errorbar((1:mk),mn,SE,'k+')
end
title([Bstr ' ordinal ISI statistics'])