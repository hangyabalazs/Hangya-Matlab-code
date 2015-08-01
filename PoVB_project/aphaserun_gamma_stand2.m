function aphaserun_gamma_stand2
%APHASERUN_GAMMA_STAND2    Phase analysis.
%   APHASERUN_GAMMA_STAND2 calculates action potential phases relative to
%   beta-band LFP oscillation. Phase histograms for all first spikes are
%   calculated. Mean angles, mean resultant length values and integrated
%   power in the 10-20 Hz band are saved in an Excel table. Group angle and
%   mean res. length statistics are saved for Po, VB, PoVPM, LD and nRT.
%   Edit program code to modify input and output directories!
%
%   APHASERUN_GAMMA_STAND2 performs EEG standarization before calculation
%   Hilbert-transform (see APHASE_GAMMA).
%
%   See also APHASE_GAMMA, ASTANORM, WATSON2 and WATSONTWOFIT2.

% Directories
global DATAPATH
inpdir1 = 'X:\In_Vivo\_analysis\acsady_csoport\auj_istvan\mat_ket_xyl\';    % mat files
inpdir2 = [DATAPATH 'Andi\Ketxyl\Cluster\mat2\'];   % burst analysis data
resdir = [DATAPATH 'Andi\Ketxyl\Spindle_stand\'];
tabledir = ['Y:\_Projects\AUJ_ISTVAN\TABLES\'];
mm = pwd;

% Filelist
[files1 files_short1] = filelist(inpdir1);
[files2 files_short2] = filelist2(inpdir2);
files_short = intersect(files_short1,files_short2);
sf = length(files_short);

% Read Excel file
fn = [tabledir 'tablazat_balazsnak.xls'];
headerrows = 0;
[mtx ntx atx] = xlsread(fn);
ntx(1:headerrows,:) = [];
atx(1:headerrows,:) = [];

% Main
ANG = struct('Po',[],'VB',[],'PoVPM',[],'LD',[],'nRT',[]);
MVL = struct('Po',[],'VB',[],'PoVPM',[],'LD',[],'nRT',[]);
for o = 1:sf
    fname = files_short{o};     % load
    cmps = strread(fname,'%s','delimiter','_.');
    fname2 = [cmps{1} '_' cmps{2}];
    ff = [inpdir1 fname];
    load(ff)
    eeg0 = data(:,2)';
    eeg = eeg0(1:20:end);    % downsample on 1000 Hz
    clear data eeg0
    sr = 1000;
    ff2 = [inpdir2 fname(1:end-4) '_CLUST2.mat'];
    load(ff2)
    
    fst = Burst(1,:);
    vb1 = vdisc(fst);
    vb1 = round(vb1/20);    % burst first spikes
    
    ssi = vdisc;       % allfirstspikes
    for k = 1:size(Burst,2)
        ssi(Burst(1,k)+1:Burst(2,k)) = 0;
    end
    afsp = ssi(ssi>0);
    afsp = round(afsp/20);    % downsample unit on 1000 Hz
    
    [aang_afsp dinx_afsp] = aphase_gamma(eeg,afsp,sr);    % PHASE - all first spikes
    n_afsp = length(aang_afsp);
    ftm_afsp = sum(exp(1).^(i*aang_afsp)) / n_afsp;    % first trigonometric moment
    ang_afsp = angle(ftm_afsp);   % mean angle
    mvl_afsp = abs(ftm_afsp);     % mean resultant length
    aang_afsp = aang_afsp * 180 / pi;
    ang_afsp = ang_afsp * 180 / pi;
    edges = [-180:20:180];
    cnts = (edges(1:end-1) + edges(2:end)) / 2;
    [nm_afsp,xout_afsp] = histc(aang_afsp,edges);   % phase histogram
    nm_afsp = nm_afsp(1:end-1);
    
    n = n_afsp;
    z = n * (mvl_afsp ^ 2);  % Rayleigh's Z statistic
    p = exp(1) ^ (-1 * z) * (1 + (2 * z - z ^ 2) / ...
        (4 * n) - (24 * z - 132 * z ^ 2 + 76 * z ^ 3 - 9 * z ^ 4) / (288 * n ^ 2));
    sl = 0.01;  % level of significance
    siglev = -1 * log(sl);
    if p > sl
        ang_afsp = NaN;
        mvl_afsp = NaN;
    end
    
    eeg2 = eeg(1:5:end);        % downsample on 200 Hz
    [y,w] = b_fft2(eeg2,200);     % FFT
    inxs = w > 10 & w < 20;
    inxs2 = w > 5 & w < 25;
    integrated_power = sum(y(inxs));
    H = figure;
    plot(w(inxs2),y(inxs2))
    hold on
    plot(w(inxs2),smooth(y(inxs2),100,'moving'),'r')
    ylim([0 500000])
    x_lim = xlim;
    y_lim = ylim;
    text((x_lim(2)-x_lim(1))/2+x_lim(1),(y_lim(2)-y_lim(1))*0.75+y_lim(1),...
        [num2str(integrated_power/1000000) ' * 10^6'])
    cd(resdir)
    saveas(H,[fname '_fft.fig'])
    close(H)
    allpower = sum(y(2:end));
    py = y / allpower;
    H = figure;
    plot(w(inxs2),py(inxs2))
    hold on
    plot(w(inxs2),smooth(py(inxs2),100,'moving'),'r')
    ylim([0 10^(-4)])
    x_lim = xlim;
    y_lim = ylim;
    text((x_lim(2)-x_lim(1))/2+x_lim(1),(y_lim(2)-y_lim(1))*0.75+y_lim(1),num2str(integrated_power/allpower))
    cd(resdir)
    saveas(H,[fname '_fftnorm.fig'])
    close(H)
    
    inx = [find(strcmp({atx{:,1}},fname2)) find(strcmp({atx{:,2}},fname2))];
    loc = atx{inx,3};     % fill new Excel table
    atx2{inx,1} = ntx{inx,1};
    atx2{inx,2} = ntx{inx,2};
    atx2{inx,3} = ntx{inx,3};
    atx2{inx,4} = ang_afsp;
    atx2{inx,5} = mvl_afsp;
    atx2{inx,6} = integrated_power / 1000000;
    atx2{inx,7} = integrated_power / allpower;
    
    eval(['ANG.' loc '(end+1) = ang_afsp;']);     % group statistics
    eval(['MVL.' loc '(end+1) = mvl_afsp;']);
end

% Save
cd(resdir)
H = mvlbox(MVL);    % group mvl statistics (Mann-Whitney U-test)
title('MVL')
saveas(H,'MVL.fig');
saveas(H,'MVL.eps');

H = angfig(ANG);    % group angle statistics (Watson-test)
title('ANG')
saveas(H,'ANG.fig');
saveas(H,'ANG.eps');

header = {'bas' 'bic' 'localization' 'angle' 'MRL' 'beta power *10^(-6)' 'normalized beta power'};
xlswrite('betatable',header,'Sheet1','A1');     % write Excel file
xlswrite('betatable',atx2,'Sheet1','A2');     % write Excel file
cd(mm)

% -------------------------------------------------------------------------
function [files2 files2_short] = filelist(inpdir)

% List of filenames
files = dir(inpdir);
files = files(3:end);
files2 = struct('name',[],'date',[],'bytes',[],'isdir',[]);
files2_short = {};
for i = 1:length(files)
    if ~files(i).isdir
        files2(end+1) = files(i);
        files2_short{end+1} = files(i).name;
    end
end
files2 = files2(2:end);

% -------------------------------------------------------------------------
function [files2 files2_short] = filelist2(inpdir)

% List of filenames
files = dir(inpdir);
files = files(3:end);
files2 = struct('name',[],'date',[],'bytes',[],'isdir',[]);
files2_short = {};
for i = 1:length(files)
    if ~files(i).isdir
        files2(end+1) = files(i);
        files2_short{end+1} = [files(i).name(1:end-11) '.mat'];
    end
end
files2 = files2(2:end);

% -------------------------------------------------------------------------
function H = mvlbox(B)

m1 = B.Po(~isnan(B.Po));
m2 = B.VB(~isnan(B.VB));
m3 = B.PoVPM(~isnan(B.PoVPM));
m4 = B.LD(~isnan(B.LD));
m5 = B.nRT(~isnan(B.nRT));

H = figure;
boxplot([m1 m2 m3 m4 m5],[zeros(size(m1)) ones(size(m2)) 2*ones(size(m3))...
    3*ones(size(m4)) 4*ones(size(m5))],'labels',[{'Po'} {'VB'} {'PoVPM'}...
    {'LD'} {'nRT'}]);
[Wp_eu,Wh_eu] = ranksum(m1,m2,'alpha',0.05);
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
function H = angfig(B)

m1 = B.Po(~isnan(B.Po));
m2 = B.VB(~isnan(B.VB));
m3 = B.PoVPM(~isnan(B.PoVPM));
m4 = B.LD(~isnan(B.LD));
m5 = B.nRT(~isnan(B.nRT));
cm1 = circular_mean(m1,'deg');
cm2 = circular_mean(m2,'deg');
cm3 = circular_mean(m3,'deg');
cm4 = circular_mean(m4,'deg');
cm5 = circular_mean(m5,'deg');
se1 = circular_SE(m1,'deg');
se2 = circular_SE(m2,'deg');
se3 = circular_SE(m3,'deg');
se4 = circular_SE(m4,'deg');
se5 = circular_SE(m5,'deg');

H = figure;
bar(1:5,[cm1 cm2 cm3 cm4 cm5])
hold on
errorbar(1:5,[cm1 cm2 cm3 cm4 cm5],[se1 se2 se3 se4 se5],'b+')
set(gca,'XTickLabel',[{'Po'} {'VB'} {'PoVPM'} {'LD'} {'nRT'}])
[U2 Wp] = b_watsontwo(m1',m2');
if Wp < 0.05
    Wh = 1;
else
    Wh = 0;
end
if Wh
    clr = 'red';
else
    clr = 'black';
end
y_lim = ylim;
x_lim = xlim;
tpos1 = x_lim(1) + (x_lim(2) - x_lim(1)) / 5;
tpos2 = y_lim(1) + (y_lim(2)-y_lim(1)) * 9 / 10;
text(tpos1,tpos2,num2str(Wp),'Color',clr,'Horizontalalignment','center')