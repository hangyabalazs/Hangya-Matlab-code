function czphase2
%CZPHASE2   Theta phase for hippocampal interneurons.
%   CZPHASE2 plots and saves theta phase histograms for hippocampal
%   interneurons. Theta phase is calculated via Hilbert-transform of the EEG
%   filtered between 4 and 12 Hz. First trigonometric moments are also
%   saved for all cells.
%
%   This version should be used for files derived from .Ntt using 
%   CZSPIKESORT!
%
%   See also CZPHASE and CZSPIKESORT.

% Directories
global DATAPATH
inpdir_eeg = [DATAPATH 'Czurko\czurko_EEG\'];
inpdir_unit = [DATAPATH 'Czurko\discriminated2\'];
resdir = [DATAPATH 'Czurko\czurko_EEG\phase\resubmission\pos\'];
mm = pwd;
cd(resdir)

% Load
xlsname = [inpdir_eeg 'EEG3c.xls'];
[ntz mtz] = xlsread(xlsname,'pos');
sf = size(mtz,1);   % number of pairs
ftms = zeros(1,sf);
ftms_rayleigh = [];
ftms_rao = [];
PRayleigh = [];
PRao = [];
for o = 1:sf
    disp(o)
    dr = [inpdir_unit mtz{o,4} '\' mtz{o,1} '\' mtz{o,5} '\'];
    fn = [dr 'clusters.mat'];
    load(fn)        % load unit
    cl = mtz{o,2};
    cmps = strread(cl,'%s','delimiter','_');
    cell_code = str2num(cmps{end});
    corrv = ntz(o);
    inx = find(CellNumbers==cell_code);
    vdisc = (TimeStamps(inx) - TimeStamps(1)) / 1000000 + corrv;
    vdisc = vdisc';
    vdisc = vdisc(vdisc>0);
    
    fn = [inpdir_eeg 'EEG_' mtz{o,3} '_' mtz{o,1} '.mat'];
    load(fn)        % load EEG
    eval(['Eeg = ' mtz{o,3} ';']);
    eval(['clear ' mtz{o,3}]);
    eeg = Eeg.values;
    sr = 1 / Eeg.interval;
    eeg_start = Eeg.start;
%     eeg_start = 0;
    eeg_end = eeg_start + (Eeg.length - 1) * Eeg.interval;
    eeg_times = eeg_start:Eeg.interval:eeg_end;
    
% Downsample EEG
    eeg2 = eeg(1:5:end);
    eeg_times2 = eeg_times(1:5:end);
    sr2 = sr / 5;
        
% Phase calculation
    [hang, hmvl, ftm, angs] = thetaphase(eeg2,vdisc,sr2,eeg_times2);
    ftms(o) = ftm;
    [Z,p_rayleigh,U,p_rao] = b_rao(angs);
    if p_rayleigh < 0.001
        ftms_rayleigh = [ftms_rayleigh ftm];
    end
    PRayleigh = [PRayleigh p_rayleigh];
    if p_rao(2) <= 0.001
        ftms_rao = [ftms_rao ftm];
    end
    PRao = [PRao; p_rao];
    
    edges = -180:20:180;     % edges for phase histogram
    cnts = (edges(1:end-1) + edges(2:end)) / 2;
    [nm,xout] = histc(angs*180/pi,edges);   % phase histogram
    nm = nm(1:end-1);
    H = figure;
    B = bar(cnts,nm'/length(angs));
    set(B,'FaceColor',[0.16 0.38 0.27])
    y_lim = ylim;
    axis([-200 200 y_lim(1) y_lim(2)])
    str = ['\it{Mean angle: }' '\bf ' num2str(hang*180/pi)];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','red')
    str = ['\it{Mean resultant length: }' '\bf ' num2str(hmvl)];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/12,str,'Color','red')
    str = ['\it{n: }' '\bf ' num2str(length(angs))];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/6,str,'Color','red')
    ts = [mtz{o,1} ' ' mtz{o,2}];
    ts(ts=='_') = ' ';
    title(ts)
    
% Save
    fns = [ts '_PHASEHIST.fig'];
    saveas(H,fns)
    fns = [ts '_STAT.mat'];
    save(fns,'p_rayleigh','p_rao')
    fns = [ts '_ANGS.mat'];
    save(fns,'angs','p_rayleigh')
end
fns = ['FTMS_pos_rayleigh.mat'];
save(fns,'ftms_rayleigh')
fns = ['FTMS_pos_rao.mat'];
save(fns,'ftms_rao')
fns = ['PVAL.mat'];
save(fns,'PRayleigh','PRao')
cd(mm)



% -------------------------------------------------------------------------
function [hang, hmvl, ftm, bahee] = thetaphase(eeg,vdisc,sr,eeg_times)

% Filtering EEG
nqf = sr / 2;
flt = fir1(4096,[4 12]/nqf,'band');      % bandpass filtering
feeg = filtfilt(flt,1,eeg);

% Hilbert transformation of the EEG
ahee = angle(hilbert(feeg));
aheedeg = ahee * (180 / pi);

% Check criteria:
% 1. discard cicles with EEG amp. lower then 2SD
% 2. discard cicles shorter then 83 ms
% 3. discard cicles longer then 250 ms
fn = find(-diff(ahee)>2*pi-0.3);
sd = std(feeg);
inx = [];
for k = 1:length(fn)-1
    seeg = feeg(fn(k):fn(k+1));
    axs = max(seeg) - min(seeg);
    sahee = ahee(fn(k):fn(k+1));
    stim = eeg_times(fn(k):fn(k+1));
    if (axs < 2 * sd)  || (fn(k+1) - fn(k) < 0.083 * sr) || (fn(k+1) - fn(k) > 0.25 * sr)
        inx = [inx; find(vdisc>stim(1)&vdisc<stim(end))];
    end
end
vdisc(inx) = [];

% Phase angles - Hilbert
vdisc(vdisc<eeg_times(1)|vdisc>eeg_times(end)) = [];
lvd = length(vdisc);
bahee = zeros(1,lvd);
for k = 1:lvd    % interpolate
    inx1 = find(eeg_times<vdisc(k),1,'last');
    inx2 = inx1 + 1;
    rto = (vdisc(k) - eeg_times(inx1)) / (eeg_times(inx2) - eeg_times(inx1));
    if rto > 1
        error('Technical error 55!')
    end
    bahee(k) = ahee(inx1) + (ahee(inx2) - ahee(inx1)) * rto;
end
n = length(bahee);
ftm = sum(exp(1).^(i*bahee)) / n;    % first trigonometric moment
hang = angle(ftm);   % mean angle
hmvl = abs(ftm);     % mean resultant length