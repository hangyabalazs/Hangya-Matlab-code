function arestrictstat_phase3
%ARESTRICTSTAT_PHASE3   Mean phase histograms for bounded frequency band.
%   ARESTRICTSTAT_PHASE3 calculates summary phase statistics for Po, VPM,
%   PoVPM, LD and nRT data. It visualizes circular mean and mean resultant
%   length for all cells on a polar plot. Edit code to modify input and
%   output directories!
%
%   See also ARESTRICTSTAT_BURST.

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
ftm_as = struct('Po',[],'VB',[],'PoVPM',[],'LD',[],'nRT',[]);
ftm_fs = struct('Po',[],'VB',[],'PoVPM',[],'LD',[],'nRT',[]);
ftm_sp = struct('Po',[],'VB',[],'PoVPM',[],'LD',[],'nRT',[]);
ftm_afsp = struct('Po',[],'VB',[],'PoVPM',[],'LD',[],'nRT',[]);
dsc = 100;      % 100 element of each phase sample
for o = 3:length(dr)
    inpadd = dr(o).name;    % load burst data
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
        eval(['[ftm_as.' loc '(end+1)] = mvlmn(aang_as/180*pi);']);
    end
    if n_fs >= dsc
        eval(['[ftm_fs.' loc '(end+1)] = mvlmn(aang_fs/180*pi);']);
    end
    if n_sp >= dsc
       eval(['[ftm_sp.' loc '(end+1)] = mvlmn(aang_sp/180*pi);']);
    end
    if n_afsp >= dsc
        eval(['[ftm_afsp.' loc '(end+1)] = mvlmn(aang_afsp/180*pi);']);
    end
    
end

% Plot and save
dbclear if error
cd(resdir)

H = figure;
subplot(2,2,1)
ftm = [0 ftm_as.Po];
CMP = compass(real(exp(1).^(i*0.1)),imag(exp(1).^(i*0.1)));
set(CMP,'Color','white')
hold on
CMP = [];
for k = 2:length(ftm)
    CMP(end+1) = compass(real(ftm(k)),imag(ftm(k)));
    set(CMP(end),'Color','black','LineWidth',2)
end
title(gca,['Po all spikes'])
subplot(2,2,2)
ftm = [0 ftm_afsp.Po];
CMP = compass(real(exp(1).^(i*0.1)),imag(exp(1).^(i*0.1)));
set(CMP,'Color','white')
hold on
CMP = [];
for k = 2:length(ftm)
    CMP(end+1) = compass(real(ftm(k)),imag(ftm(k)));
    set(CMP(end),'Color','black','LineWidth',2)
end
ftm_Po = abs(ftm(2:end));
title(gca,['Po all first spikes'])
subplot(2,2,3)
ftm = [0 ftm_fs.Po];
CMP = compass(real(exp(1).^(i*0.1)),imag(exp(1).^(i*0.1)));
set(CMP,'Color','white')
hold on
CMP = [];
for k = 2:length(ftm)
    CMP(end+1) = compass(real(ftm(k)),imag(ftm(k)));
    set(CMP(end),'Color','black','LineWidth',2)
end
title(gca,['Po burst first spikes'])
subplot(2,2,4)
ftm = [0 ftm_sp.Po];
CMP = compass(real(exp(1).^(i*0.1)),imag(exp(1).^(i*0.1)));
set(CMP,'Color','white')
hold on
CMP = [];
for k = 2:length(ftm)
    CMP(end+1) = compass(real(ftm(k)),imag(ftm(k)));
    set(CMP(end),'Color','black','LineWidth',2)
end
title(gca,['Po single spikes'])
saveas(H,'compass_Po.fig')
saveas(H,'compass_Po.eps')

H = figure;
subplot(2,2,1)
ftm = [0 ftm_as.VB];
CMP = compass(real(exp(1).^(i*0.1)),imag(exp(1).^(i*0.1)));
set(CMP,'Color','white')
hold on
CMP = [];
for k = 2:length(ftm)
    CMP(end+1) = compass(real(ftm(k)),imag(ftm(k)));
    set(CMP(end),'Color','black','LineWidth',2)
end
title(gca,['VB all spikes'])
subplot(2,2,2)
ftm = [0 ftm_afsp.VB];
CMP = compass(real(exp(1).^(i*0.1)),imag(exp(1).^(i*0.1)));
set(CMP,'Color','white')
hold on
CMP = [];
for k = 2:length(ftm)
    CMP(end+1) = compass(real(ftm(k)),imag(ftm(k)));
    set(CMP(end),'Color','black','LineWidth',2)
end
ftm_VB = abs(ftm(2:end));
title(gca,['VB all first spikes'])
subplot(2,2,3)
ftm = [0 ftm_fs.VB];
CMP = compass(real(exp(1).^(i*0.1)),imag(exp(1).^(i*0.1)));
set(CMP,'Color','white')
hold on
CMP = [];
for k = 2:length(ftm)
    CMP(end+1) = compass(real(ftm(k)),imag(ftm(k)));
    set(CMP(end),'Color','black','LineWidth',2)
end
title(gca,['VB burst first spikes'])
subplot(2,2,4)
ftm = [0 ftm_sp.VB];
CMP = compass(real(exp(1).^(i*0.1)),imag(exp(1).^(i*0.1)));
set(CMP,'Color','white')
hold on
CMP = [];
for k = 2:length(ftm)
    CMP(end+1) = compass(real(ftm(k)),imag(ftm(k)));
    set(CMP(end),'Color','black','LineWidth',2)
end
title(gca,['VB single spikes'])
saveas(H,'compass_VB.fig')
saveas(H,'compass_VB.eps')

H = figure;
subplot(2,2,1)
ftm = [0 ftm_as.PoVPM];
CMP = compass(real(exp(1).^(i*0.1)),imag(exp(1).^(i*0.1)));
set(CMP,'Color','white')
hold on
CMP = [];
for k = 2:length(ftm)
    CMP(end+1) = compass(real(ftm(k)),imag(ftm(k)));
    set(CMP(end),'Color','black','LineWidth',2)
end
title(gca,['PoVPM all spikes'])
subplot(2,2,2)
ftm = [0 ftm_afsp.PoVPM];
CMP = compass(real(exp(1).^(i*0.1)),imag(exp(1).^(i*0.1)));
set(CMP,'Color','white')
hold on
CMP = [];
for k = 2:length(ftm)
    CMP(end+1) = compass(real(ftm(k)),imag(ftm(k)));
    set(CMP(end),'Color','black','LineWidth',2)
end
ftm_PoVPM = abs(ftm(2:end));
title(gca,['PoVPM all first spikes'])
subplot(2,2,3)
ftm = [0 ftm_fs.PoVPM];
CMP = compass(real(exp(1).^(i*0.1)),imag(exp(1).^(i*0.1)));
set(CMP,'Color','white')
hold on
CMP = [];
for k = 2:length(ftm)
    CMP(end+1) = compass(real(ftm(k)),imag(ftm(k)));
    set(CMP(end),'Color','black','LineWidth',2)
end
title(gca,['PoVPM burst first spikes'])
subplot(2,2,4)
ftm = [0 ftm_sp.PoVPM];
CMP = compass(real(exp(1).^(i*0.1)),imag(exp(1).^(i*0.1)));
set(CMP,'Color','white')
hold on
CMP = [];
for k = 2:length(ftm)
    CMP(end+1) = compass(real(ftm(k)),imag(ftm(k)));
    set(CMP(end),'Color','black','LineWidth',2)
end
title(gca,['PoVPM single spikes'])
saveas(H,'compass_PoVPM.fig')
saveas(H,'compass_PoVPM.eps')

H = figure;
subplot(2,2,1)
ftm = [0 ftm_as.LD];
CMP = compass(real(exp(1).^(i*0.1)),imag(exp(1).^(i*0.1)));
set(CMP,'Color','white')
hold on
CMP = [];
for k = 2:length(ftm)
    CMP(end+1) = compass(real(ftm(k)),imag(ftm(k)));
    set(CMP(end),'Color','black','LineWidth',2)
end
title(gca,['LD all spikes'])
subplot(2,2,2)
ftm = [0 ftm_afsp.LD];
CMP = compass(real(exp(1).^(i*0.1)),imag(exp(1).^(i*0.1)));
set(CMP,'Color','white')
hold on
CMP = [];
for k = 2:length(ftm)
    CMP(end+1) = compass(real(ftm(k)),imag(ftm(k)));
    set(CMP(end),'Color','black','LineWidth',2)
end
ftm_LD = abs(ftm(2:end));
title(gca,['LD all first spikes'])
subplot(2,2,3)
ftm = [0 ftm_fs.LD];
CMP = compass(real(exp(1).^(i*0.1)),imag(exp(1).^(i*0.1)));
set(CMP,'Color','white')
hold on
CMP = [];
for k = 2:length(ftm)
    CMP(end+1) = compass(real(ftm(k)),imag(ftm(k)));
    set(CMP(end),'Color','black','LineWidth',2)
end
title(gca,['LD burst first spikes'])
subplot(2,2,4)
ftm = [0 ftm_sp.LD];
CMP = compass(real(exp(1).^(i*0.1)),imag(exp(1).^(i*0.1)));
set(CMP,'Color','white')
hold on
CMP = [];
for k = 2:length(ftm)
    CMP(end+1) = compass(real(ftm(k)),imag(ftm(k)));
    set(CMP(end),'Color','black','LineWidth',2)
end
title(gca,['LD single spikes'])
saveas(H,'compass_LD.fig')
saveas(H,'compass_LD.eps')

H = figure;
subplot(2,2,1)
ftm = [0 ftm_as.nRT];
CMP = compass(real(exp(1).^(i*0.1)),imag(exp(1).^(i*0.1)));
set(CMP,'Color','white')
hold on
CMP = [];
for k = 2:length(ftm)
    CMP(end+1) = compass(real(ftm(k)),imag(ftm(k)));
    set(CMP(end),'Color','black','LineWidth',2)
end
title(gca,['nRT all spikes'])
subplot(2,2,2)
ftm = [0 ftm_afsp.nRT];
CMP = compass(real(exp(1).^(i*0.1)),imag(exp(1).^(i*0.1)));
set(CMP,'Color','white')
hold on
CMP = [];
for k = 2:length(ftm)
    CMP(end+1) = compass(real(ftm(k)),imag(ftm(k)));
    set(CMP(end),'Color','black','LineWidth',2)
end
ftm_nRT = abs(ftm(2:end));
title(gca,['nRT all first spikes'])
subplot(2,2,3)
ftm = [0 ftm_fs.nRT];
CMP = compass(real(exp(1).^(i*0.1)),imag(exp(1).^(i*0.1)));
set(CMP,'Color','white')
hold on
CMP = [];
for k = 2:length(ftm)
    CMP(end+1) = compass(real(ftm(k)),imag(ftm(k)));
    set(CMP(end),'Color','black','LineWidth',2)
end
title(gca,['nRT burst first spikes'])
subplot(2,2,4)
ftm = [0 ftm_sp.nRT];
CMP = compass(real(exp(1).^(i*0.1)),imag(exp(1).^(i*0.1)));
set(CMP,'Color','white')
hold on
CMP = [];
for k = 2:length(ftm)
    CMP(end+1) = compass(real(ftm(k)),imag(ftm(k)));
    set(CMP(end),'Color','black','LineWidth',2)
end
title(gca,['nRT single spikes'])
saveas(H,'compass_nRT.fig')
saveas(H,'compass_nRT.eps')

H = figure;     % MVL-s
plot(ftm_Po,'.')
title('Po')
saveas(H,'mvls_Po.fig')
saveas(H,'mvls_Po.eps')
plot(ftm_VB,'.')
title('VB')
saveas(H,'mvls_VB.fig')
saveas(H,'mvls_VB.eps')
plot(ftm_PoVPM,'.')
title('PoVPM')
saveas(H,'mvls_PoVPM.fig')
saveas(H,'mvls_PoVPM.eps')
plot(ftm_LD,'.')
title('LD')
saveas(H,'mvls_LD.fig')
saveas(H,'mvls_LD.eps')
plot(ftm_nRT,'.')
title('nRT')
saveas(H,'mvls_nRT.fig')
saveas(H,'mvls_nRT.eps')

H = figure;
plot(sort(ftm_Po))
title('Po')
saveas(H,'mvlsort_Po.fig')
saveas(H,'mvlsort_Po.eps')
plot(sort(ftm_VB))
title('VB')
saveas(H,'mvlsort_VB.fig')
saveas(H,'mvlsort_VB.eps')
plot(sort(ftm_PoVPM))
title('PoVPM')
saveas(H,'mvlsort_PoVPM.fig')
saveas(H,'mvlsort_PoVPM.eps')
plot(sort(ftm_LD))
title('LD')
saveas(H,'mvlsort_LD.fig')
saveas(H,'mvlsort_LD.eps')
plot(sort(ftm_nRT))
title('nRT')
saveas(H,'mvlsort_nRT.fig')
saveas(H,'mvlsort_nRT.eps')

mvl_Po = mean(ftm_Po);
mvl_VB = mean(ftm_VB);
mvl_PoVPM = mean(ftm_PoVPM);
mvl_LD = mean(ftm_LD);
mvl_nRT = mean(ftm_nRT);
save mvl_means mvl_Po mvl_VB mvl_PoVPM mvl_LD mvl_nRT

cd(mm)

% -------------------------------------------------------------------------
function ftm = mvlmn(fr1_all)

ftm = sum(exp(1).^(i*fr1_all)) / length(fr1_all);    % first trigonometric moment