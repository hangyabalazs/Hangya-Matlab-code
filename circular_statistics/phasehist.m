function [hang hmvl ftm bahee ahee inx H Hc] = phasehist(eeg,vdisc,sr,fr1,fr2,fo,mlt,ip)
%PHASEHIST   Phase histogram.
%   [A R FTM ANGS AHEE INX] = PHASEHIST(EEG,VD,SR,FR1,FR2,FO,MLT,IP)
%   calculates phase histograms for neuronal spikes (VD) relative to EEG
%   sampled at SR. Phase is calculated via Hilbert-transform of the EEG
%   filtered between FR1 and FR2 Hz and standardized. Filtering is
%   performed by FIR filter applying bidirectional zero-phase filtering
%   with a filter order of FO (FO, optional input parameter; default, 2048).
%   Output arguments are mean phase (A), mean vector length (R), first
%   trigonometric moment (FTM), all unit phase values (ANGS), all EEG phase
%   values (AHEE) and indeces of spikes left out from the calculations (see
%   criteria below) (INX). Cycles with a length out of the frequency range
%   or with an amplitude lower than mean + MLT * SD are discarded (MLT,
%   optional input parameter; default value, 2). IP, optional input
%   argument; if IP equals to 1, control plots are displayed (default, 0).
%   Control plots include filtered trace and phase series as well as phase
%   histogram for all EEG phase values. Significance of phase locking is
%   tested with Rayleigh's test and Rao's spacing test (see RAO for
%   DETAILS).
%
%   [A R FTM ANGS AHEE INX H HC] = PHASEHIST(EEG,VD,SRFR1,FR2) returns
%   figure handles of unit phase histogram and a control phase histogram,
%   from all EEG phase values.
%
%   See also RAO, SOMPHASECALL4, APHASE_STAND and CZPHASE2.

%   Balazs Hangya
%   Cold Spring Harbor Laboratory
%   1 Bungtown Road, Cold Spring Harbor
%   New York 11724, USA
%   balazs.cshl@gmail.com
%
%   Institute of Experimental Medicine
%   Szigony street 43, Budapest 1083, Hungary
%   hangyab@koki.hu

% Input argument check
error(nargchk(5,8,nargin))
if nargin < 8
    ip = 0;
end
if nargin < 7
    mlt = 2;
end
if nargin < 6
    fo = 2048;
end

% Filtering EEG
nqf = sr / 2;
if fr1 == 0
    flt = fir1(fo,fr2/nqf,'low');
elseif fr2 == Inf
    flt = fir1(fo,fr1/nqf,'high');
else
    flt = fir1(fo,[fr1 fr2]/nqf,'band');      % bandpass filtering
end
feeg = filtfilt(flt,1,eeg);
feeg = standardize(feeg);

% Hilbert transformation of the EEG
ahee = angle(hilbert(feeg));
aheedeg = ahee * (180 / pi);
if ip
    figure
    plot(standardize(eeg))
    hold on
    plot(feeg,'g')
    plot(ahee,'c')
end

% Check criteria:
% 1. discard cicles with EEG amp. lower then MLT * SD
% 2. discard cicles shorter then 1/fr2 s
% 3. discard cicles longer then 1/fr1 s
fn0 = valuecrossing(1:length(ahee),ahee(:)',0,'down');
fn = round(fn0);
sd = std(feeg);
inx = [];
vdisc(vdisc<0|vdisc>length(eeg)) = [];
for k = 1:length(fn)-1
    seeg = feeg(fn(k)+1:fn(k+1));
    axs = max(seeg) - min(seeg);
    if (axs < mlt * sd)  || (fn(k+1) - fn(k) > 1 / fr1 * sr) || (fn(k+1) - fn(k) < 1 / fr2 * sr)
        inx = [inx find(vdisc>fn0(k)&vdisc<fn0(k+1))]; %#ok<AGROW>
        if ip
            plot(fn(k)+1:fn(k+1),seeg,'r')
        end
    end
end
vdisc(inx) = [];

% Phase angles - Hilbert
bahee = ahee(round(vdisc));
n = length(bahee);
ftm = sum(exp(1).^(1i*bahee)) / n;    % first trigonometric moment
hang = angle(ftm);   % mean angle
hmvl = abs(ftm);     % mean resultant length

% Plot phase histogram
angs = bahee;
edges = -180:20:180;     % edges for phase histogram
cnts = (edges(1:end-1) + edges(2:end)) / 2;
nm = histc(angs*180/pi,edges);   % phase histogram
nm = nm(1:end-1);
H = figure;
B = bar(cnts,nm'/length(angs));
set(B,'FaceColor',[0.16 0.38 0.27])
y_lim = ylim;
axis([-200 200 y_lim(1) y_lim(2)])
str = ['\it{Mean angle: }' '\bf ' num2str(hang*180/pi)];
text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','black')
str = ['\it{Mean resultant length: }' '\bf ' num2str(hmvl)];
text(60,y_lim(2)-1.5*(y_lim(2)-y_lim(1))/6,str,'Color','black')
str = ['\it{n: }' '\bf ' num2str(length(angs))];
text(60,y_lim(2)-2*(y_lim(2)-y_lim(1))/6,str,'Color','black')
if ~isempty(bahee)      % significance testing
    [Z,p_rayleigh,U,p_rao] = b_rao(bahee(:)');
    if p_rayleigh < 0.001
        clr_ray = 'red';
    else
        clr_ray = 'black';
    end
    if p_rao(2) <= 0.05
        clr_rao = 'red';
    else
        clr_rao = 'black';
    end
    str = ['\it{Raylegh''s p = }' '\bf ' num2str(p_rayleigh)];
    text(60,y_lim(2)-2.5*(y_lim(2)-y_lim(1))/6,str,'Color',clr_ray)
    str = ['\it{Rao''s p < }' '\bf ' num2str(p_rao(2))];
    text(60,y_lim(2)-3*(y_lim(2)-y_lim(1))/6,str,'Color',clr_rao)
end

% Plot phase histogram of all EEG phase values
if ip || nargout > 7
    nm = histc(ahee*180/pi,edges);   % control phase histogram (all EEG phase values)
    nm = nm(1:end-1);
    Hc = figure;
    B = bar(cnts,nm'/length(ahee));
    set(B,'FaceColor',[0.16 0.38 0.27])
    y_lim = ylim;
    axis([-200 200 y_lim(1) y_lim(2)])
end