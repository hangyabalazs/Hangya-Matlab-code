function somphasecall4
%SOMPHASECALL4   Phase calculations.
%   SOMPHASECALL4 adopts the Klausberger-type ripple detection algorithm to
%   beta episodes.
%
%   See also SOMPHASECALL2, SOMPHASECALL3 and SSPW.

% Directories
global DATAPATH
inpdir = '\\science\Kepecs\Hyun\For Balazs\';
resdir = [DATAPATH '\SOM\Phase5\'];
mm = pwd;
cd(resdir)

% Phase histogram
dr = dir(inpdir);
for k = 4:length(dr)
    
    if (dr(k).isdir && strncmp(dr(k).name,'.',1)==0)
        
        % Load EEG
        inpadd = dr(k).name;
        drr = [inpdir inpadd];
        [fl0 fl] = filelist(drr);
        fst = find(strncmp(fl,'CSC',3));
        eegname = fl{fst(1)};
        eegnm = eegname(1:end-4);
        load([inpdir inpadd filesep eegname])
        load([inpdir inpadd filesep 'EVENTS'])  %load events
        eval(['frs = ' eegnm '_SampleFrequencies;'])    % sampling rate
        eval(['eegts = ' eegnm '_TimeStamps;'])    % EEG time stamps
        eval(['clear ' eegnm '_TimeStamps'])    % EEG time stamps
        if max(frs) > 1020
            continue
        end
        
        % Remove LaserStim part
        eval(['eeg0 = ' eegnm '_Samples;']);  % EEG
        eval(['clear ' eegnm '_Samples']);  % EEG
        sr = 1010;
        dt = 1 / sr;
        eeg0 = eeg0(:);
        time = (0:length(eeg0)-1) * dt + eegts(1)/1000000;     % EEG time stamps
        t0 = eegts(1)/1000000;
        if Events_TimeStamps(1) > 1e8
            Events_TimeStamps = Events_TimeStamps / 1000000;
        end
        idx = find(strcmp(Events_EventStrings, 'Cheetah 32 DIO board 0 input TTL (0x0200)'),1);
        idx_start = find(time>Events_TimeStamps(idx),1,'first');
        eeg = eeg0(1:idx_start-1);
        clear eeg0 Events_TimeStamps
        
        % Resample
        eeg = resample(eeg,99,100);
        sr = 1000;
        dt = 1 / sr;
%         time = (0:length(eeg)-1) * dt + eegts(1) / 1000000;     % EEG time stamps
        
        
        %         %plot spectrogram
        %         spectrogram1(eeg,dt);
        %         fname = ['Spectrogram_' inpadd '_' eegnm '.fig'];
        %         saveas(gcf,fname)
        
        % Load unit
        fst0 = find(strncmp(fl,'TT',2));
        for cells = 2:2
            unitname = fl{fst0(cells)};
            unitnm = unitname(1:end-4);
            load([inpdir inpadd filesep unitname])
            TS = TS - eegts(1)/1000000;
            TS = TS(TS<(length(eeg)/sr));
            vdisc = TS' * sr;
            
            % Phase calculation
            [hang hmvl ftm bahee ahee Hhist Hctrl] = somphase(eeg,vdisc,sr,3,10);   % theta
            fnm = ['THETA_PHASE_' inpadd '_' eegnm '_' unitnm '.fig'];
            saveas(Hhist,fnm)
%             fnm = ['THETA_PHASEctrl_' inpadd '_' eegnm '_' unitnm '.fig'];
%             saveas(Hctrl,fnm)
            if ~isempty(bahee)
                [Z,p_rayleigh,U,p_rao] = b_rao(bahee');
            end
            disp(p_rao)
            %fnm = ['THETA_PHASE_' inpadd '_' eegnm '_' unitnm '.mat'];
            %save(fnm,'hang','hmvl','ftm','bahee')
            
            betaL = 20;     % beta band boundaries
            betaU = 35;
            [bg fin feeg] = sspw(eeg,sr,betaL,betaU,3);   % thresholding on filtered signal
            lenb = length(bg);
            nbinx = [];
            Heeg = figure;
            plot(eeg)
            hold on
            plot(feeg,'r')
            for bts = 1:lenb    % restrict spike train to detected segments
                ind1 = bg(bts);
                ind2 = fin(bts);
                nbinx = [nbinx find(vdisc>ind1&vdisc<ind2)];
                plot(ind1:ind2,feeg(ind1:ind2),'g')
            end
            vdisc2 = vdisc(nbinx);
            feeg = (feeg - mean(feeg)) / std(feeg);     % standardize
            ahee = angle(hilbert(feeg));
            bahee = ahee(round(vdisc2));
            
            n = length(bahee);      % cirular statistics
            ftm = sum(exp(1).^(i*bahee)) / n;    % first trigonometric moment
            hang = angle(ftm);   % mean angle
            hmvl = abs(ftm);     % mean resultant length
            
            angs = bahee;       % plot phase histogram
            edges = -180:20:180;     % edges for phase histogram
            cnts = (edges(1:end-1) + edges(2:end)) / 2;
            [nm,xout] = histc(angs*180/pi,edges);   % phase histogram
            nm = nm(1:end-1);
            Hhist = figure;
            B = bar(cnts,nm'/length(angs));
            set(B,'FaceColor',[0.16 0.38 0.27])
            y_lim = ylim;
            axis([-200 200 y_lim(1) y_lim(2)])
            str = ['\it{Mean angle: }' '\bf ' num2str(hang*180/pi)];
            text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','red')
            str = ['\it{Mean resultant length: }' '\bf ' num2str(hmvl)];
            text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/12,str,'Color','red')
            str = ['\it{n: }' '\bf ' num2str(length(angs))];
            text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/6,str,'Color','red')%cd('beta')
            
            fnm = ['BETA_PHASE_' inpadd '_' eegnm '_' unitnm '.fig'];
            saveas(Hhist,fnm)
            %fnm = ['BETA_PHASE_' inpadd '_' eegnm '_' unitnm '.mat'];
            %save(fnm,'hang','hmvl','ftm','bahee')
            %cd ..
            clear ahee
            
            close all
        end
    end
end
cd(mm)



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
function [hang hmvl ftm bahee ahee H Hc] = somphase(eeg,vdisc,sr,fr1,fr2)
%SOMPHASE   Phase histogram.
%   [A R FTM ANGS AHEE] = SOMPHASE(EEG,VD,SRFR1,FR2) calculates phase
%   histograms for neuronal spikes (VD) relative to EEG sampled at SR.
%   Phase is calculated via Hilbert-transform of the EEG filtered between
%   FR1 and FR2 Hz. Output arguments are mean phase (A), mean vector length
%   (R), first trigonometric moment (FTM), all unit phase values (ANGS) and
%   all EEG phase values (AHEE).
%
%   [A R FTM ANGS AHEE H HC] = SOMPHASE(EEG,VD,SRFR1,FR2) returns figure
%   handles of unit phase histogram and a control phase histogram, from all
%   EEG phase values.
%
%   See also CZPHASE2.

% Filtering EEG
nqf = sr / 2;
flt = fir1(2048,[fr1 fr2]/nqf,'band');      % bandpass filtering
feeg = filtfilt(flt,1,eeg);
feeg = (feeg - mean(feeg)) / std(feeg);
figure
plot((eeg-mean(eeg))/std(eeg))
hold on
plot(feeg,'g')

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
        plot(fn(k)+1:fn(k+1),seeg,'r')
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

[nm,xout] = histc(ahee*180/pi,edges);   % control phase histogram (all EEG phase values)
nm = nm(1:end-1);
Hc = figure;
B = bar(cnts,nm'/length(ahee));
set(B,'FaceColor',[0.16 0.38 0.27])
y_lim = ylim;
axis([-200 200 y_lim(1) y_lim(2)])

% -------------------------------------------------------------------------
function [st nd feeg] = sspw(eeg,sr,betaL,betaU,mlt)
%SSPW   Beta detection.
%   [ST ND FEEG] = SSPW(EEG,SR,L,U,M) detectects beta epochs (between L and
%   U Hz) in EEG sampled at SR and returns starting and end points of
%   epochs in ST and ND, respectively. SSPW uses the following criteria for
%   beta epochs: root-mean-square of filtered EEG should reach mean(RMS) +
%   std(RMS) and peak RMS should reach mean(RMS) + M * std(RMS). Filtered
%   EEG is returned in FEEG.
%
%   See also SSPW_NONTHETA.

% Input argument check
error(nargchk(5,5,nargin))

% Filtering
nqf = sr / 2;
b = fir1(2048,[betaL betaU]/nqf);
feeg = filtfilt(b,1,eeg);

% Root mean square
leneeg = length(feeg);
wl1 = 100;
lle1 = floor(leneeg/wl1) * wl1;
feeg21 = reshape(feeg(1:lle1),wl1,lle1/wl1);
rms1 = sqrt(sum(feeg21.^2)) / sqrt(wl1);
wl2 = 100;
lle2 = floor(leneeg/wl2) * wl2;
feeg22 = reshape(feeg(1:lle2),wl2,lle2/wl2);
rms2 = sqrt(sum(feeg22.^2)) / sqrt(wl2);

% Discriminate RMS: RMS peak during sharpwave should reach mean(RMS) + mlt * std(RMS)
mrms1 = mean(rms1);
sdrms1 = std(rms1);
thr = mrms1 + mlt * sdrms1;
pks = disc(rms1,thr);

% Ripple start, end: RMS should cross mean(RMS) + std(RMS)
lenp = length(pks);
st = zeros(1,lenp);
nd = zeros(1,lenp);
mrms2 = mean(rms2);
sdrms2 = std(rms2);
v = mrms2 + sdrms2;
lup = rms2 < v & [rms2(2:end) v] > v;   % value-crossings
lup2 = find(lup) * wl2;
ldown = rms2 < v & [v rms2(1:end-1)] > v;
ldown2 = find(ldown) * wl2;
pks = pks * wl1;
for k = 1:length(pks)
    sst = lup2(find(lup2<pks(k),1,'last'));
    nnd = ldown2(find(ldown2>pks(k),1,'first'));
    if isempty(sst) || isempty(nnd)
        continue
    end
    st(k) = sst;
    nd(k) = nnd;
end
st = st(st>0);
nd = nd(nd>0);
[st m1] = unique(st);
[nd m2] = unique(nd);
if ~isequal(m1,m2)
    error('Technical error 50.')
end