function [hang, hmvl, ftm, bahee] = somphase(eeg,vdisc,sr,eeg_times,fr1,fr2)
%SOMPHASE   Phase histogram.
%   [A R FTM ANGS] = SOMPHASE(EEG,VD,SR,EEG_TIMES,FR1,FR2) calculates phase
%   histograms for neuronal spikes (VD) relative to EEG sampled at SR,
%   given at time instances stored in EEG_TIMES. Phase is calculated via
%   Hilbert-transform of the EEG filtered between FR1 and FR2 Hz. Output
%   arguments are mean phase (A), mean vector length (R), first
%   trigonometric moment (FTM) and all phase values (ANGS).
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
sahees = [];
for k = 1:length(fn)-1
    seeg = feeg(fn(k)+1:fn(k+1));
    axs = max(seeg) - min(seeg);
    sahee = ahee(fn(k)+1:fn(k+1));
    stim = linterp(1:length(eeg_times),eeg_times,[fn0(k) fn0(k+1)]);
    if (axs < 2 * sd)  || (fn(k+1) - fn(k) > 1 / fr1 * sr) || (fn(k+1) - fn(k) < 1 / fr2 * sr)
        inx = [inx; find(vdisc>stim(1)&vdisc<stim(2))];
    else
        sahees = [sahees; sahee];
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
    if ahee(inx1) - ahee(inx2) < pi
        bahee(k) = ahee(inx1) + (ahee(inx2) - ahee(inx1)) * rto;
    else
        bahee(k) = ahee(inx1) - 2 * pi + (ahee(inx2) - (ahee(inx1) - 2 * pi)) * rto;
    end
%     if abs(vdisc(k)-eeg_times(inx1)) < abs(vdisc(k)-eeg_times(inx2))
%         bahee(k) = ahee(inx1);
%     else
%         bahee(k) = ahee(inx2);
%     end
end
bahee(bahee<-pi) = bahee(bahee<-pi) + 2 * pi;
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