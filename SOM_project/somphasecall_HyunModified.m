function somphasecall4
%SOMPHASECALL2   Phase calculations.
%   SOMPHASECALL2 is functionally similar to SOMPHASECALL but uses an EEG
%   resampling approach instead of interpolation, which results in faster
%   execution.
%
%   See also SOMPHASECALL.

% Directories
global DATAPATH
inpdir = '\\science\Kepecs\Hyun\For Balazs\';
resdir = [DATAPATH '\SOM\Phase3\'];

mm = pwd;
cd(resdir)

% Phase histogram
dr = dir(inpdir);
for k = 1:length(dr)
    
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
        if max(frs) > 1020
            continue
        end
        
        % Remove LaserStim part
        eval(['eeg0 = ' eegnm '_Samples;']);  % EEG
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
        
        % Resample
        eeg = resample(eeg,99,100);
        sr = 1000;
        dt = 1 / sr;
        time = (0:length(eeg)-1) * dt + eegts(1) / 1000000;     % EEG time stamps
        
        
        %         %plot spectrogram
        %         spectrogram1(eeg,dt);
        %         fname = ['Spectrogram_' inpadd '_' eegnm '.fig'];
        %         saveas(gcf,fname)
        
        % Load unit
        fst0 = find(strncmp(fl,'TT',2));
        for cells = 1:length(fst0)
            unitname = fl{fst0(cells)};
            unitnm = unitname(1:end-4);
            load([inpdir inpadd filesep unitname])
            TS = TS - eegts(1)/1000000;
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
            
%             [hang, hmvl, ftm, bahee] = somphase(eeg,TS'*sr,sr,20,35);   % beta
%             %cd('beta')
%             fnm = ['BETA_PHASE_' inpadd '_' eegnm '_' unitnm '.fig'];
%             saveas(gcf,fnm)
%             %fnm = ['BETA_PHASE_' inpadd '_' eegnm '_' unitnm '.mat'];
%             %save(fnm,'hang','hmvl','ftm','bahee')
%             %cd ..
            
            % Spike Triggered Average
            wn = 2 * sr;    % 2 sec. window
            [StaTheta StaIndexTheta1 StaIndexTheta2 nn] = astanorm(round(vdisc),eeg,wn);
            titlestr = [inpadd ' ' eegnm ' ' unitnm];
            titlestr(titlestr=='_') = ' ';
            H = stafig(StaTheta,StaIndexTheta1,StaIndexTheta2,nn,wn,sr,titlestr);
            fnm = ['THETA_STA_' inpadd '_' eegnm '_' unitnm '.fig'];
            saveas(H,fnm);
            
            % Autocorrelation
            wn = sr;
            [H1 H2 trsc] = somccg(TS,TS,wn);  % CCG
            fnm = ['THETA_AUTO_' inpadd '_' eegnm '_' unitnm '.fig'];
            saveas(H1,fnm);
            
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
% figure
% plot((eeg-mean(eeg))/std(eeg))
% hold on
% plot(feeg,'g')

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
%         plot(fn(k)+1:fn(k+1),seeg,'r')
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
function H = stafig(sta,sta_index1,sta_index2,nn,wn,sr,titlestr)

time = linspace(-wn/sr/2,wn/sr/2,length(sta));
H = figure;
plot(time,sta,'LineWidth',1.5)
ach = allchild(H);     % figure title
ax = findobj(ach,'type','axes');
title(ax(end),titlestr)
x_lim = xlim;
y_lim = ylim;
str = ['\it{Max-mean: }' '\bf ' num2str(sta_index1)];
text(0.1,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','red')
str = ['\it{Max: }' '\bf ' num2str(sta_index2)];
text(0.1,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/12,str,'Color','red')
str = ['\it{n: }' '\bf ' num2str(nn)];
text(0.1,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/6,str,'Color','red')