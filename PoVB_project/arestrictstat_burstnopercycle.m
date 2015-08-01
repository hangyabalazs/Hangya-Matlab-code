function arestrictstat_burstnopercycle
%ARESTRICTSTAT_BURSTNOPERCYCLE   Burst number per cycle histograms for bounded frequency band.
%   ARESTRICTSTAT_BURSTNOPERCYCLE calculates summary phase statistics for
%   Po, VPM, PoVPM, LD and nRT data. Original burst no. per cycle
%   histograms are downsampled (to 100 elements) and pooled to provide
%   balanced sets for tests. Summary histograms and box plots of percentage
%   of missed cycles are generated ans saved. Edit code to modify input and
%   output directories!
%
%   See also ARESTRICTSTAT_BURST and AFRE_RESTRICT.

% Input argument check
error(nargchk(0,0,nargin))
dbstop if error

% Directories
global DATAPATH
inpdir = 'Y:\_Projects\AUJ_ISTVAN\DATA\MAT\mat_ket_xyl\';
inpdir2 = [DATAPATH 'Andi\Ketxyl\FreBandRestrict_phase_stand\'];   % phase analysis data
tabledir = ['Y:\_Projects\AUJ_ISTVAN\TABLES\'];
resdir = [DATAPATH 'Andi\Ketxyl\RestrictStat_stand\'];
mm = pwd;
dr = dir(inpdir);

% Main
BNOPC = struct('Po',{{}},'VB',{{}},'PoVPM',{{}},'LD',{{}},'nRT',{{}});
CM = struct('Po',[],'VB',[],'PoVPM',[],'LD',[],'nRT',[]);
dsc = 100;      % 100 element of each sample
for o = 3:length(dr)
    inpadd = dr(o).name;    % load phase data
    cd(inpdir)
    cd(inpadd)
    cd('bas')
    ddr = dir(pwd);
    fn = ddr(end).name(1:end-4);
    cmps = strread(fn,'%s','delimiter','_');
    fname = [cmps{1} '_' cmps{2}];
    ff = [inpdir2 fn '_PHASE.mat'];
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
    
    n_cycnb = length(cycnb);     % burst number per cycle
    if n_cycnb >= dsc
        rp = randperm(n_cycnb);
        cycnb_ds = cycnb(rp(1:dsc));
        eval(['BNOPC.' loc ' = [BNOPC.' loc ' cycnb_ds];']);
        cycmis = length(find(cycnb==0)) / length(cycnb);
        eval(['CM.' loc '(end+1) = cycmis;']);
    end
end

% Plot and save
num.Po = length(cell2mat(BNOPC.Po)) / dsc;
num.VB = length(cell2mat(BNOPC.VB)) / dsc;
num.PoVPM = length(cell2mat(BNOPC.PoVPM)) / dsc;
num.LD = length(cell2mat(BNOPC.LD)) / dsc;
num.nRT = length(cell2mat(BNOPC.nRT)) / dsc;
BNOPCpo = cell2mat(BNOPC.Po);
BNOPCvb = cell2mat(BNOPC.VB);
BNOPCpovpm = cell2mat(BNOPC.PoVPM);
BNOPCld = cell2mat(BNOPC.LD);
BNOPCnrt = cell2mat(BNOPC.nRT);
[nm_po xout_po] = hist(BNOPCpo(BNOPCpo>0),(1:10));     % distr. of burst no./cycle
[nm_vb xout_vb] = hist(BNOPCvb(BNOPCvb>0),(1:10));
[nm_povpm xout_povpm] = hist(BNOPCpovpm(BNOPCpovpm>0),(1:10));
[nm_ld xout_ld] = hist(BNOPCld(BNOPCld>0),(1:10));
[nm_nrt xout_nrt] = hist(BNOPCnrt(BNOPCnrt>0),(1:10));

dbclear if error
cd(resdir)

H = figure;
B = bar(xout_po,nm_po);
set(B,'FaceColor',[0.16 0.38 0.27])
title(gca,['Po distr. of burst no./cycle'])
saveas(H,'burstnopercycle_Po.fig')
saveas(H,'burstnopercycle_Po.eps')

H = figure;
B = bar(xout_vb,nm_vb);
set(B,'FaceColor',[0.16 0.38 0.27])
title(gca,['VB distr. of burst no./cycle'])
saveas(H,'burstnopercycle_VB.fig')
saveas(H,'burstnopercycle_VB.eps')

H = figure;
B = bar(xout_povpm,nm_povpm);
set(B,'FaceColor',[0.16 0.38 0.27])
title(gca,['PoVPM distr. of burst no./cycle'])
saveas(H,'burstnopercycle_PoVPM.fig')
saveas(H,'burstnopercycle_PoVPM.eps')

H = figure;
B = bar(xout_ld,nm_ld);
set(B,'FaceColor',[0.16 0.38 0.27])
title(gca,['LD distr. of burst no./cycle'])
saveas(H,'burstnopercycle_LD.fig')
saveas(H,'burstnopercycle_LD.eps')

H = figure;
B = bar(xout_nrt,nm_nrt);
set(B,'FaceColor',[0.16 0.38 0.27])
title(gca,['nRT distr. of burst no./cycle'])
saveas(H,'burstnopercycle_nRT.fig')
saveas(H,'burstnopercycle_nRT.eps')

H = burstbox(CM);
title('Cycles missed')
saveas(H,'Missedcycles_boxplot.fig');
saveas(H,'Missedcycles_boxplot.eps');

cd(mm)

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