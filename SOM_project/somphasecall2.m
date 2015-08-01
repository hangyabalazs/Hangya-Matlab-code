function somphasecall2
%SOMPHASECALL2   Phase calculations.
%   SOMPHASECALL2 is functionally similar to SOMPHASECALL but uses an EEG
%   resampling approach instead of interpolation, which results in faster
%   execution.
%
%   See also SOMPHASECALL.

% Directories
global DATAPATH
inpdir = '\\science\Kepecs\recordings\for Balazs\new\';
resdir = [DATAPATH 'SOM\Phase2\'];
mm = pwd;
cd(resdir)

% Phase histogram
dr = dir(inpdir);
for k = 1:length(dr)
    if dr(k).isdir
        
        % Load EEG
        inpadd = dr(k).name;
        drr = [inpdir inpadd];
        [fl0 fl] = filelist(drr);
        fst = find(strncmp(fl,'CSC',3));
        eegname = fl{fst(1)};
        eegnm = eegname(1:end-4);
        load([inpdir inpadd '\' eegname])
        eval(['frs = ' eegnm '_SampleFrequencies;'])    % sampling rate
        eval(['eegts = ' eegnm '_TimeStamps;'])    % EEG time stamps
        if max(frs) > 1020
            continue
        end
        eval(['eeg = ' eegnm '_Samples;'])  % EEG
        eeg = eeg(:);
        eeg = resample(eeg,99,100);
        sr = 1000;
        dt = 1 / sr;
        time = (0:length(eeg)-1)*dt+eegts(1)/1000000;     % EEG time stamps
        
        % Load unit
        fst = find(strncmp(fl,'TT',2));
        for cells = 1:length(fst)
            unitname = fl{fst(cells)};
            unitnm = unitname(1:end-4);
            load([inpdir inpadd '\' unitname])
            TS = TS - eegts(1)/1000000;
            
            % Beta alignment
            
            % Phase calculation
            [hang, hmvl, ftm, bahee] = somphase(eeg,TS'*sr,sr,3,10);   % theta
            cd('theta')
            fnm = ['THETA_PHASE_' inpadd '_' eegnm '_' unitnm '.fig'];
            saveas(gcf,fnm)
            fnm = ['THETA_PHASE_' inpadd '_' eegnm '_' unitnm '.mat'];
            save(fnm,'hang','hmvl','ftm','bahee')
            cd ..
            
            [hang, hmvl, ftm, bahee] = somphase(eeg,TS'*sr,sr,20,35);   % beta
            cd('beta')
            fnm = ['BETA_PHASE_' inpadd '_' eegnm '_' unitnm '.fig'];
            saveas(gcf,fnm)
            fnm = ['BETA_PHASE_' inpadd '_' eegnm '_' unitnm '.mat'];
            save(fnm,'hang','hmvl','ftm','bahee')
            cd ..
            close all
        end
    end
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
        files2_short{end+1} = files(i).name;
    end
end
files2 = files2(2:end);

% -------------------------------------------------------------------------
function [hang, hmvl, ftm, bahee] = somphase(eeg,vdisc,sr,fr1,fr2)
%SOMPHASE   Phase histogram.
%   [A R FTM ANGS] = SOMPHASE(EEG,VD,SRFR1,FR2) calculates phase histograms
%   for neuronal spikes (VD) relative to EEG sampled at SR. Phase is
%   calculated via Hilbert-transform of the EEG filtered between FR1 and
%   FR2 Hz. Output arguments are mean phase (A), mean vector length (R),
%   first trigonometric moment (FTM) and all phase values (ANGS).
%
%   See also CZPHASE2.

% Filtering EEG
nqf = sr / 2;
flt = fir1(4096,[fr1 fr2]/nqf,'band');      % bandpass filtering
feeg = filtfilt(flt,1,eeg);
feeg = (feeg - mean(feeg)) / std(feeg);

% Hilbert transformation of the EEG
ahee = angle(hilbert(feeg));
aheedeg = ahee * (180 / pi);

% Check criteria:
% 1. discard cicles with EEG amp. lower then 2SD
% 2. discard cicles shorter then 1/fr2 s
% 3. discard cicles longer then 1/fr1 s
% fn = find(-diff(ahee)>2*pi-0.3);
fn0 = valuecrossing(1:length(ahee),ahee',0,'down');
fn = round(fn0);
sd = std(feeg);
inx = [];
vdisc(vdisc<0|vdisc>length(eeg)) = [];
for k = 1:length(fn)-1
    seeg = feeg(fn(k)+1:fn(k+1));
    axs = max(seeg) - min(seeg);
    if (axs < 2 * sd)  || (fn(k+1) - fn(k) > 1 / fr1 * sr) || (fn(k+1) - fn(k) < 1 / fr2 * sr)
        inx = [inx find(vdisc>fn0(k)&vdisc<fn0(k+1))];
    end
end
vdisc(inx) = [];

% Phase angles - Hilbert
bahee = ahee(round(vdisc));
n = length(bahee);
ftm = sum(exp(1).^(i*bahee)) / n;    % first trigonometric moment
hang = angle(ftm);   % mean angle
hmvl = abs(ftm);     % mean resultant length

% Plot phase histogram
angs = bahee;
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