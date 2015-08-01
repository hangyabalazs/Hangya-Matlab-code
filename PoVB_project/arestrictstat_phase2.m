function arestrictstat_phase2
%ARESTRICTSTAT_PHASE2   Mean phase histograms for bounded frequency band.
%   ARESTRICTSTAT_PHASE2 calculates summary phase statistics for Po, VPM,
%   PoVPM, LD and nRT data. Original phase histograms are downsampled (to
%   100 elements) and pooled to provide balanced sets for circular tests.
%   Edit code to modify input and output directories!
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
resdir = [DATAPATH 'Andi\Ketxyl\RestrictStat_stand_newtest\'];
mm = pwd;
dr = dir(inpdir);

% Main
AS = struct('Po',{{}},'VB',{{}},'PoVPM',{{}},'LD',{{}},'nRT',{{}});
FS = struct('Po',{{}},'VB',{{}},'PoVPM',{{}},'LD',{{}},'nRT',{{}});
SP = struct('Po',{{}},'VB',{{}},'PoVPM',{{}},'LD',{{}},'nRT',{{}});
AFSP = struct('Po',{{}},'VB',{{}},'PoVPM',{{}},'LD',{{}},'nRT',{{}});
edges = -180:20:180;     % edges for phase histogram
cnts = (edges(1:end-1) + edges(2:end)) / 2;
dsc = 100;      % 100 element of each phase sample
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
    
    n_fs = length(aang_fs);     % burst first spikes
    n_sp = length(aang_sp);     % single spikes
    n_as = length(aang_as);     % all spikes
    n_afsp = length(aang_afsp);     % all first spikes
    if n_as >= dsc
        rp = randperm(n_as);
        aang_as_ds = aang_as(rp(1:dsc));
        eval(['AS.' loc ' = [AS.' loc ' aang_as_ds];']);
    end
    if n_fs >= dsc
        rp = randperm(n_fs);
        aang_fs_ds = aang_fs(rp(1:dsc));
        eval(['FS.' loc ' = [FS.' loc ' aang_fs_ds];']);
    end
    if n_sp >= dsc
        rp = randperm(n_sp);
        aang_sp_ds = aang_sp(rp(1:dsc));
        eval(['SP.' loc ' = [SP.' loc ' aang_sp_ds];']);
    end
    if n_afsp >= dsc
        rp = randperm(n_afsp);
        aang_afsp_ds = aang_afsp(rp(1:dsc));
        eval(['AFSP.' loc ' = [AFSP.' loc ' aang_afsp_ds];']);
    end
end

% Plot and save
num_as.Po = length(cell2mat(AS.Po)) / dsc;
num_as.VB = length(cell2mat(AS.VB)) / dsc;
num_as.PoVPM = length(cell2mat(AS.PoVPM)) / dsc;
num_as.LD = length(cell2mat(AS.LD)) / dsc;
num_as.nRT = length(cell2mat(AS.nRT)) / dsc;
nm_as_po = histc(cell2mat(AS.Po),edges) / (num_as.Po * dsc);   % phase histogram
nm_as_po = nm_as_po(1:end-1);
alla_as.Po = nm_as_po;
nm_as_vb = histc(cell2mat(AS.VB),edges) / (num_as.VB * dsc);
nm_as_vb = nm_as_vb(1:end-1);
alla_as.VB = nm_as_vb;
nm_as_povpm = histc(cell2mat(AS.PoVPM),edges) / (num_as.PoVPM * dsc);
nm_as_povpm = nm_as_povpm(1:end-1);
alla_as.PoVPM = nm_as_povpm;
nm_as_ld = histc(cell2mat(AS.LD),edges) / (num_as.LD * dsc);
nm_as_ld = nm_as_ld(1:end-1);
alla_as.LD = nm_as_ld;
nm_as_nrt = histc(cell2mat(AS.nRT),edges) / (num_as.nRT * dsc);
nm_as_nrt = nm_as_nrt(1:end-1);
alla_as.nRT = nm_as_nrt;

num_fs.Po = length(cell2mat(FS.Po)) / dsc;
num_fs.VB = length(cell2mat(FS.VB)) / dsc;
num_fs.PoVPM = length(cell2mat(FS.PoVPM)) / dsc;
num_fs.LD = length(cell2mat(FS.LD)) / dsc;
num_fs.nRT = length(cell2mat(FS.nRT)) / dsc;
nm_fs_po = histc(cell2mat(FS.Po),edges) / (num_fs.Po * dsc);   % phase histogram
nm_fs_po = nm_fs_po(1:end-1);
alla_fs.Po = nm_fs_po;
nm_fs_vb = histc(cell2mat(FS.VB),edges) / (num_fs.VB * dsc);
nm_fs_vb = nm_fs_vb(1:end-1);
alla_fs.VB = nm_fs_vb;
nm_fs_povpm = histc(cell2mat(FS.PoVPM),edges) / (num_fs.PoVPM * dsc);
nm_fs_povpm = nm_fs_povpm(1:end-1);
alla_fs.PoVPM = nm_fs_povpm;
nm_fs_ld = histc(cell2mat(FS.LD),edges) / (num_fs.LD * dsc);
nm_fs_ld = nm_fs_ld(1:end-1);
alla_fs.LD = nm_fs_ld;
nm_fs_nrt = histc(cell2mat(FS.nRT),edges) / (num_fs.nRT * dsc);
nm_fs_nrt = nm_fs_nrt(1:end-1);
alla_fs.nRT = nm_fs_nrt;

num_sp.Po = length(cell2mat(SP.Po)) / dsc;
num_sp.VB = length(cell2mat(SP.VB)) / dsc;
num_sp.PoVPM = length(cell2mat(SP.PoVPM)) / dsc;
num_sp.LD = length(cell2mat(SP.LD)) / dsc;
num_sp.nRT = length(cell2mat(SP.nRT)) / dsc;
nm_sp_po = histc(cell2mat(SP.Po),edges) / (num_sp.Po * dsc);   % phase histogram
nm_sp_po = nm_sp_po(1:end-1);
alla_sp.Po = nm_sp_po;
nm_sp_vb = histc(cell2mat(SP.VB),edges) / (num_sp.VB * dsc);
nm_sp_vb = nm_sp_vb(1:end-1);
alla_sp.VB = nm_sp_vb;
nm_sp_povpm = histc(cell2mat(SP.PoVPM),edges) / (num_sp.PoVPM * dsc);
nm_sp_povpm = nm_sp_povpm(1:end-1);
alla_sp.PoVPM = nm_sp_povpm;
nm_sp_ld = histc(cell2mat(SP.LD),edges) / (num_sp.LD * dsc);
nm_sp_ld = nm_sp_ld(1:end-1);
alla_sp.LD = nm_sp_ld;
nm_sp_nrt = histc(cell2mat(SP.nRT),edges) / (num_sp.nRT * dsc);
nm_sp_nrt = nm_sp_nrt(1:end-1);
alla_sp.nRT = nm_sp_nrt;

num_afsp.Po = length(cell2mat(AFSP.Po)) / dsc;
num_afsp.VB = length(cell2mat(AFSP.VB)) / dsc;
num_afsp.PoVPM = length(cell2mat(AFSP.PoVPM)) / dsc;
num_afsp.LD = length(cell2mat(AFSP.LD)) / dsc;
num_afsp.nRT = length(cell2mat(AFSP.nRT)) / dsc;
nm_afsp_po = histc(cell2mat(AFSP.Po),edges) / (num_afsp.Po * dsc);   % phase histogram
nm_afsp_po = nm_afsp_po(1:end-1);
alla_afsp.Po = nm_afsp_po;
nm_afsp_vb = histc(cell2mat(AFSP.VB),edges) / (num_afsp.VB * dsc);
nm_afsp_vb = nm_afsp_vb(1:end-1);
alla_afsp.VB = nm_afsp_vb;
nm_afsp_povpm = histc(cell2mat(AFSP.PoVPM),edges) / (num_afsp.PoVPM * dsc);
nm_afsp_povpm = nm_afsp_povpm(1:end-1);
alla_afsp.PoVPM = nm_afsp_povpm;
nm_afsp_ld = histc(cell2mat(AFSP.LD),edges) / (num_afsp.LD * dsc);
nm_afsp_ld = nm_afsp_ld(1:end-1);
alla_afsp.LD = nm_afsp_ld;
nm_afsp_nrt = histc(cell2mat(AFSP.nRT),edges) / (num_afsp.nRT * dsc);
nm_afsp_nrt = nm_afsp_nrt(1:end-1);
alla_afsp.nRT = nm_afsp_nrt;

dbclear if error
cd(resdir)

H = figure;
subplot(2,2,1)
B = bar(cnts,alla_as.Po);
set(B,'FaceColor',[0.16 0.38 0.27])
title(gca,['Po all spikes'])
y_lim = ylim;
axis([-200 200 y_lim(1) y_lim(2)])
str = ['Mean of ' num2str(num_as.Po) ' histograms.'];
text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','red')
subplot(2,2,2)
B = bar(cnts,alla_afsp.Po);
set(B,'FaceColor',[0.16 0.38 0.27])
title(gca,['Po all first spikes'])
y_lim = ylim;
axis([-200 200 y_lim(1) y_lim(2)])
str = ['Mean of ' num2str(num_afsp.Po) ' histograms.'];
text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','red')
subplot(2,2,3)
B = bar(cnts,alla_fs.Po);
set(B,'FaceColor',[0.16 0.38 0.27])
title(gca,['Po burst first spikes'])
y_lim = ylim;
axis([-200 200 y_lim(1) y_lim(2)])
str = ['Mean of ' num2str(num_fs.Po) ' histograms.'];
text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','red')
subplot(2,2,4)
B = bar(cnts,alla_sp.Po);
set(B,'FaceColor',[0.16 0.38 0.27])
title(gca,['Po single spikes'])
y_lim = ylim;
axis([-200 200 y_lim(1) y_lim(2)])
str = ['Mean of ' num2str(num_sp.Po) ' histograms.'];
text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','red')
saveas(H,'phase_Po2.fig')
saveas(H,'phase_Po2.eps')

H = figure;
subplot(2,2,1)
B = bar(cnts,alla_as.VB);
set(B,'FaceColor',[0.16 0.38 0.27])
title(gca,['VB all spikes'])
y_lim = ylim;
axis([-200 200 y_lim(1) y_lim(2)])
str = ['Mean of ' num2str(num_as.VB) ' histograms.'];
text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','red')
subplot(2,2,2)
B = bar(cnts,alla_afsp.VB);
set(B,'FaceColor',[0.16 0.38 0.27])
title(gca,['VB all first spikes'])
y_lim = ylim;
axis([-200 200 y_lim(1) y_lim(2)])
str = ['Mean of ' num2str(num_afsp.VB) ' histograms.'];
text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','red')
subplot(2,2,3)
B = bar(cnts,alla_fs.VB);
set(B,'FaceColor',[0.16 0.38 0.27])
title(gca,['VB burst first spikes'])
y_lim = ylim;
axis([-200 200 y_lim(1) y_lim(2)])
str = ['Mean of ' num2str(num_fs.VB) ' histograms.'];
text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','red')
subplot(2,2,4)
B = bar(cnts,alla_sp.VB);
set(B,'FaceColor',[0.16 0.38 0.27])
title(gca,['VB single spikes'])
y_lim = ylim;
axis([-200 200 y_lim(1) y_lim(2)])
str = ['Mean of ' num2str(num_sp.VB) ' histograms.'];
text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','red')
saveas(H,'phase_VB2.fig')
saveas(H,'phase_VB2.eps')

H = figure;
subplot(2,2,1)
B = bar(cnts,alla_as.PoVPM);
set(B,'FaceColor',[0.16 0.38 0.27])
title(gca,['PoVPM all spikes'])
y_lim = ylim;
axis([-200 200 y_lim(1) y_lim(2)])
str = ['Mean of ' num2str(num_as.PoVPM) ' histograms.'];
text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','red')
subplot(2,2,2)
B = bar(cnts,alla_afsp.PoVPM);
set(B,'FaceColor',[0.16 0.38 0.27])
title(gca,['PoVPM all first spikes'])
y_lim = ylim;
axis([-200 200 y_lim(1) y_lim(2)])
str = ['Mean of ' num2str(num_afsp.PoVPM) ' histograms.'];
text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','red')
subplot(2,2,3)
B = bar(cnts,alla_fs.PoVPM);
set(B,'FaceColor',[0.16 0.38 0.27])
title(gca,['PoVPM burst first spikes'])
y_lim = ylim;
axis([-200 200 y_lim(1) y_lim(2)])
str = ['Mean of ' num2str(num_fs.PoVPM) ' histograms.'];
text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','red')
subplot(2,2,4)
B = bar(cnts,alla_sp.PoVPM);
set(B,'FaceColor',[0.16 0.38 0.27])
title(gca,['PoVPM single spikes'])
y_lim = ylim;
axis([-200 200 y_lim(1) y_lim(2)])
str = ['Mean of ' num2str(num_sp.PoVPM) ' histograms.'];
text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','red')
saveas(H,'phase_PoVPM2.fig')
saveas(H,'phase_PoVPM2.eps')

H = figure;
subplot(2,2,1)
B = bar(cnts,alla_as.LD);
set(B,'FaceColor',[0.16 0.38 0.27])
title(gca,['LD all spikes'])
y_lim = ylim;
axis([-200 200 y_lim(1) y_lim(2)])
str = ['Mean of ' num2str(num_as.LD) ' histograms.'];
text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','red')
subplot(2,2,2)
B = bar(cnts,alla_afsp.LD);
set(B,'FaceColor',[0.16 0.38 0.27])
title(gca,['LD all first spikes'])
y_lim = ylim;
axis([-200 200 y_lim(1) y_lim(2)])
str = ['Mean of ' num2str(num_afsp.LD) ' histograms.'];
text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','red')
subplot(2,2,3)
B = bar(cnts,alla_fs.LD);
set(B,'FaceColor',[0.16 0.38 0.27])
title(gca,['LD burst first spikes'])
y_lim = ylim;
axis([-200 200 y_lim(1) y_lim(2)])
str = ['Mean of ' num2str(num_fs.LD) ' histograms.'];
text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','red')
subplot(2,2,4)
B = bar(cnts,alla_sp.LD);
set(B,'FaceColor',[0.16 0.38 0.27])
title(gca,['LD single spikes'])
y_lim = ylim;
axis([-200 200 y_lim(1) y_lim(2)])
str = ['Mean of ' num2str(num_sp.LD) ' histograms.'];
text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','red')
saveas(H,'phase_LD2.fig')
saveas(H,'phase_LD2.eps')

H = figure;
subplot(2,2,1)
B = bar(cnts,alla_as.nRT);
set(B,'FaceColor',[0.16 0.38 0.27])
title(gca,['nRT all spikes'])
y_lim = ylim;
axis([-200 200 y_lim(1) y_lim(2)])
str = ['Mean of ' num2str(num_as.nRT) ' histograms.'];
text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','red')
subplot(2,2,2)
B = bar(cnts,alla_afsp.nRT);
set(B,'FaceColor',[0.16 0.38 0.27])
title(gca,['nRT all first spikes'])
y_lim = ylim;
axis([-200 200 y_lim(1) y_lim(2)])
str = ['Mean of ' num2str(num_afsp.nRT) ' histograms.'];
text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','red')
subplot(2,2,3)
B = bar(cnts,alla_fs.nRT);
set(B,'FaceColor',[0.16 0.38 0.27])
title(gca,['nRT burst first spikes'])
y_lim = ylim;
axis([-200 200 y_lim(1) y_lim(2)])
str = ['Mean of ' num2str(num_fs.nRT) ' histograms.'];
text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','red')
subplot(2,2,4)
B = bar(cnts,alla_sp.nRT);
set(B,'FaceColor',[0.16 0.38 0.27])
title(gca,['nRT single spikes'])
y_lim = ylim;
axis([-200 200 y_lim(1) y_lim(2)])
str = ['Mean of ' num2str(num_sp.nRT) ' histograms.'];
text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','red')
saveas(H,'phase_nRT2.fig')
saveas(H,'phase_nRT2.eps')

% Watson-test
phasedist1 = cell2mat(AS.Po)' / 180 * pi;
phasedist2 = cell2mat(AS.VB)' / 180 * pi;
[S.as.p S.as.mn1 S.as.mn2 S.as.mvl1 S.as.mvl2] = wtest(phasedist1,phasedist2);
phasedist1 = cell2mat(FS.Po)' / 180 * pi;
phasedist2 = cell2mat(FS.VB)' / 180 * pi;
[S.fs.p S.fs.mn1 S.fs.mn2 S.fs.mvl1 S.fs.mvl2] = wtest(phasedist1,phasedist2);
phasedist1 = cell2mat(SP.Po)' / 180 * pi;
phasedist2 = cell2mat(SP.VB)' / 180 * pi;
[S.sp.p S.sp.mn1 S.sp.mn2 S.sp.mvl1 S.sp.mvl2] = wtest(phasedist1,phasedist2);
phasedist1 = cell2mat(AFSP.Po)' / 180 * pi;
phasedist2 = cell2mat(AFSP.VB)' / 180 * pi;
[S.afsp.p S.afsp.mn1 S.afsp.mn2 S.afsp.mvl1 S.afsp.mvl2] = wtest(phasedist1,phasedist2);
phasestat_Po_VB = S;

phasedist1 = cell2mat(AS.Po)' / 180 * pi;
phasedist2 = cell2mat(AS.PoVPM)' / 180 * pi;
[S.as.p S.as.mn1 S.as.mn2 S.as.mvl1 S.as.mvl2] = wtest(phasedist1,phasedist2);
phasedist1 = cell2mat(FS.Po)' / 180 * pi;
phasedist2 = cell2mat(FS.PoVPM)' / 180 * pi;
[S.fs.p S.fs.mn1 S.fs.mn2 S.fs.mvl1 S.fs.mvl2] = wtest(phasedist1,phasedist2);
phasedist1 = cell2mat(SP.Po)' / 180 * pi;
phasedist2 = cell2mat(SP.PoVPM)' / 180 * pi;
[S.sp.p S.sp.mn1 S.sp.mn2 S.sp.mvl1 S.sp.mvl2] = wtest(phasedist1,phasedist2);
phasedist1 = cell2mat(AFSP.Po)' / 180 * pi;
phasedist2 = cell2mat(AFSP.PoVPM)' / 180 * pi;
[S.afsp.p S.afsp.mn1 S.afsp.mn2 S.afsp.mvl1 S.afsp.mvl2] = wtest(phasedist1,phasedist2);
phasestat_Po_PoVPM = S;

phasedist1 = cell2mat(AS.Po)' / 180 * pi;
phasedist2 = cell2mat(AS.LD)' / 180 * pi;
[S.as.p S.as.mn1 S.as.mn2 S.as.mvl1 S.as.mvl2] = wtest(phasedist1,phasedist2);
phasedist1 = cell2mat(FS.Po)' / 180 * pi;
phasedist2 = cell2mat(FS.LD)' / 180 * pi;
[S.fs.p S.fs.mn1 S.fs.mn2 S.fs.mvl1 S.fs.mvl2] = wtest(phasedist1,phasedist2);
phasedist1 = cell2mat(SP.Po)' / 180 * pi;
phasedist2 = cell2mat(SP.LD)' / 180 * pi;
[S.sp.p S.sp.mn1 S.sp.mn2 S.sp.mvl1 S.sp.mvl2] = wtest(phasedist1,phasedist2);
phasedist1 = cell2mat(AFSP.Po)' / 180 * pi;
phasedist2 = cell2mat(AFSP.LD)' / 180 * pi;
[S.afsp.p S.afsp.mn1 S.afsp.mn2 S.afsp.mvl1 S.afsp.mvl2] = wtest(phasedist1,phasedist2);
phasestat_Po_LD = S;

phasedist1 = cell2mat(AS.VB)' / 180 * pi;
phasedist2 = cell2mat(AS.PoVPM)' / 180 * pi;
[S.as.p S.as.mn1 S.as.mn2 S.as.mvl1 S.as.mvl2] = wtest(phasedist1,phasedist2);
phasedist1 = cell2mat(FS.VB)' / 180 * pi;
phasedist2 = cell2mat(FS.PoVPM)' / 180 * pi;
[S.fs.p S.fs.mn1 S.fs.mn2 S.fs.mvl1 S.fs.mvl2] = wtest(phasedist1,phasedist2);
phasedist1 = cell2mat(SP.VB)' / 180 * pi;
phasedist2 = cell2mat(SP.PoVPM)' / 180 * pi;
[S.sp.p S.sp.mn1 S.sp.mn2 S.sp.mvl1 S.sp.mvl2] = wtest(phasedist1,phasedist2);
phasedist1 = cell2mat(AFSP.VB)' / 180 * pi;
phasedist2 = cell2mat(AFSP.PoVPM)' / 180 * pi;
[S.afsp.p S.afsp.mn1 S.afsp.mn2 S.afsp.mvl1 S.afsp.mvl2] = wtest(phasedist1,phasedist2);
phasestat_VB_PoVPM = S;

phasedist1 = cell2mat(AS.VB)' / 180 * pi;
phasedist2 = cell2mat(AS.LD)' / 180 * pi;
[S.as.p S.as.mn1 S.as.mn2 S.as.mvl1 S.as.mvl2] = wtest(phasedist1,phasedist2);
phasedist1 = cell2mat(FS.VB)' / 180 * pi;
phasedist2 = cell2mat(FS.LD)' / 180 * pi;
[S.fs.p S.fs.mn1 S.fs.mn2 S.fs.mvl1 S.fs.mvl2] = wtest(phasedist1,phasedist2);
phasedist1 = cell2mat(SP.VB)' / 180 * pi;
phasedist2 = cell2mat(SP.LD)' / 180 * pi;
[S.sp.p S.sp.mn1 S.sp.mn2 S.sp.mvl1 S.sp.mvl2] = wtest(phasedist1,phasedist2);
phasedist1 = cell2mat(AFSP.VB)' / 180 * pi;
phasedist2 = cell2mat(AFSP.LD)' / 180 * pi;
[S.afsp.p S.afsp.mn1 S.afsp.mn2 S.afsp.mvl1 S.afsp.mvl2] = wtest(phasedist1,phasedist2);
phasestat_VB_LD = S;

phasedist1 = cell2mat(AS.PoVPM)' / 180 * pi;
phasedist2 = cell2mat(AS.LD)' / 180 * pi;
[S.as.p S.as.mn1 S.as.mn2 S.as.mvl1 S.as.mvl2] = wtest(phasedist1,phasedist2);
phasedist1 = cell2mat(FS.PoVPM)' / 180 * pi;
phasedist2 = cell2mat(FS.LD)' / 180 * pi;
[S.fs.p S.fs.mn1 S.fs.mn2 S.fs.mvl1 S.fs.mvl2] = wtest(phasedist1,phasedist2);
phasedist1 = cell2mat(SP.PoVPM)' / 180 * pi;
phasedist2 = cell2mat(SP.LD)' / 180 * pi;
[S.sp.p S.sp.mn1 S.sp.mn2 S.sp.mvl1 S.sp.mvl2] = wtest(phasedist1,phasedist2);
phasedist1 = cell2mat(AFSP.PoVPM)' / 180 * pi;
phasedist2 = cell2mat(AFSP.LD)' / 180 * pi;
[S.afsp.p S.afsp.mn1 S.afsp.mn2 S.afsp.mvl1 S.afsp.mvl2] = wtest(phasedist1,phasedist2);
phasestat_PoVPM_LD = S;
save phasestat phasestat_Po_VB phasestat_Po_PoVPM phasestat_Po_LD...
    phasestat_VB_PoVPM phasestat_VB_LD phasestat_PoVPM_LD


cd(mm)

% -------------------------------------------------------------------------
function [p mn1 mn2 mvl1 mvl2] = wtest(fr1_all,fr2_all)

[p H] = kappacompare2(fr1_all',fr2_all','rad');
ftm = sum(exp(1).^(i*fr1_all)) / length(fr1_all);    % first trigonometric moment
mn1 = angle(ftm);   % mean angle
mn1 = mod(mn1,2*pi) * 180 / pi;
mvl1 = abs(ftm);     % mean resultant length
ftm = sum(exp(1).^(i*fr2_all)) / length(fr2_all);    % first trigonometric moment
mn2 = angle(ftm);   % mean angle
mn2 = mod(mn2,2*pi) * 180 / pi;
mvl2 = abs(ftm);