function arestrictstat_phase
%ARESTRICTSTAT_PHASE   Mean phase histograms for bounded frequency band.
%   ARESTRICTSTAT_PHASE calculates summary phase statistics for Po, VPM,
%   PoVPM, LD and nRT data. Edit code to modify input
%   and output directories!

% Input argument check
error(nargchk(0,0,nargin))
dbstop if error

% Directories
global DATAPATH
inpdir = 'Y:\_Projects\AUJ_ISTVAN\DATA\MAT\mat_ket_xyl\';
inpdir2 = [DATAPATH 'Andi\Ketxyl\FreBandRestrict_phase2\'];   % phase analysis data
tabledir = ['Y:\_Projects\AUJ_ISTVAN\TABLES\'];
resdir = [DATAPATH 'Andi\Ketxyl\RestrictStat\'];
mm = pwd;
dr = dir(inpdir);

% Main
AS = struct('Po',[],'VB',[],'PoVPM',[],'LD',[],'nRT',[]);
FS = struct('Po',[],'VB',[],'PoVPM',[],'LD',[],'nRT',[]);
SP = struct('Po',[],'VB',[],'PoVPM',[],'LD',[],'nRT',[]);
AFSP = struct('Po',[],'VB',[],'PoVPM',[],'LD',[],'nRT',[]);
edges = -180:20:180;     % edges for phase histogram
cnts = (edges(1:end-1) + edges(2:end)) / 2;
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
    [nm_fs,xout_fs] = histc(aang_fs,edges);   % phase histogram
    nm_fs = nm_fs(1:end-1);
    
    n_sp = length(aang_sp);     % single spikes
    [nm_sp,xout_sp] = histc(aang_sp,edges);   % phase histogram
    nm_sp = nm_sp(1:end-1);

    n_as = length(aang_as);     % all spikes
    [nm_as,xout_as] = histc(aang_as,edges);   % phase histogram
    nm_as = nm_as(1:end-1);

    n_afsp = length(aang_afsp);     % all first spikes
    [nm_afsp,xout_afsp] = histc(aang_afsp,edges);   % phase histogram
    nm_afsp = nm_afsp(1:end-1);
    
    if n_as >= 50
        eval(['AS.' loc '(end+1,:) = nm_as/n_as;']);
    end
    if n_fs >= 50
        eval(['FS.' loc '(end+1,:) = nm_fs/n_fs;']);
    end
    if n_sp >= 50
        eval(['SP.' loc '(end+1,:) = nm_sp/n_sp;']);
    end
    if n_afsp >= 50
        eval(['AFSP.' loc '(end+1,:) = nm_afsp/n_afsp;']);
    end
    
end

% Plot and save
alla_as.Po = mean(AS.Po,1);
alla_as.VB = mean(AS.VB,1);
alla_as.PoVPM = mean(AS.PoVPM,1);
alla_as.LD = mean(AS.LD,1);
alla_as.nRT = mean(AS.nRT,1);
num_as.Po = size(AS.Po,1);
num_as.VB = size(AS.VB,1);
num_as.PoVPM = size(AS.PoVPM,1);
num_as.LD = size(AS.LD,1);
num_as.nRT = size(AS.nRT,1);

alla_fs.Po = mean(FS.Po,1);
alla_fs.VB = mean(FS.VB,1);
alla_fs.PoVPM = mean(FS.PoVPM,1);
alla_fs.LD = mean(FS.LD,1);
alla_fs.nRT = mean(FS.nRT,1);
num_fs.Po = size(FS.Po,1);
num_fs.VB = size(FS.VB,1);
num_fs.PoVPM = size(FS.PoVPM,1);
num_fs.LD = size(FS.LD,1);
num_fs.nRT = size(FS.nRT,1);

alla_sp.Po = mean(SP.Po,1);
alla_sp.VB = mean(SP.VB,1);
alla_sp.PoVPM = mean(SP.PoVPM,1);
alla_sp.LD = mean(SP.LD,1);
alla_sp.nRT = mean(SP.nRT,1);
num_sp.Po = size(SP.Po,1);
num_sp.VB = size(SP.VB,1);
num_sp.PoVPM = size(SP.PoVPM,1);
num_sp.LD = size(SP.LD,1);
num_sp.nRT = size(SP.nRT,1);

alla_afsp.Po = mean(AFSP.Po,1);
alla_afsp.VB = mean(AFSP.VB,1);
alla_afsp.PoVPM = mean(AFSP.PoVPM,1);
alla_afsp.LD = mean(AFSP.LD,1);
alla_afsp.nRT = mean(AFSP.nRT,1);
num_afsp.Po = size(AFSP.Po,1);
num_afsp.VB = size(AFSP.VB,1);
num_afsp.PoVPM = size(AFSP.PoVPM,1);
num_afsp.LD = size(AFSP.LD,1);
num_afsp.nRT = size(AFSP.nRT,1);

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
saveas(H,'phase_Po1.fig')
saveas(H,'phase_Po1.eps')

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
saveas(H,'phase_VB1.fig')
saveas(H,'phase_VB1.eps')

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
saveas(H,'phase_PoVPM1.fig')
saveas(H,'phase_PoVPM1.eps')

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
saveas(H,'phase_LD1.fig')
saveas(H,'phase_LD1.eps')

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
saveas(H,'phase_nRT1.fig')
saveas(H,'phase_nRT1.eps')

cd(mm)