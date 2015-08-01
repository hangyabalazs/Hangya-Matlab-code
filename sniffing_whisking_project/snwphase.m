function snwphase

load('C:\Balazs\_data\HSW\R3\2008-07-09_14-46-47_R3\ALIGNED_SIGNALS.mat')

% Linearize data matrices 
PELLET_H1 = -PELLET_H1';
PELLET_SN = PELLET_SN';
PELLET_W1 = -PELLET_W1';
eeg = PELLET_H1(:);
sniff = PELLET_SN(:);
whisk = PELLET_W1(:);

% Resample
sr = 1 / (ALIGNED_TIME(2) - ALIGNED_TIME(1));
[p q] = rat(1000/sr);
eeg = resample(eeg,p,q);    % resample at 1000 Hz
sniff = resample(sniff,p,q);
whisk = resample(whisk,p,q);
sr = 1000;

% Filter
nqf = sr / 2;
flt = fir1(1024,30/nqf,'low');
fsniff = filtfilt(flt,1,sniff);     % lowpass filter sniffing at 30 Hz
flt = fir1(4096,10/nqf,'high');
fwhisk = filtfilt(flt,1,whisk);     % highpass filter whisking at 10 Hz
flt = fir1(1024,[4 12]/nqf,'band');     % filter EEG in the theta band
feeg = filtfilt(flt,1,eeg);
% figure
% plot(whisk)
% hold on
% plot(fwhisk,'r')
figure
plot(eeg)
hold on
plot(feeg,'r')

% Respiration phase
resp_phase = angle(hilbert(standardize(fsniff)));    % Hilbert-transform
% figure
% plot((sniff-mean(sniff))/std(sniff))
% hold on
% plot((fsniff-mean(fsniff))/std(fsniff),'r')
% plot(resp_phase,'c')
figure
plot((sniff-mean(sniff))/std(sniff))
hold on
plot(resp_phase,'r')

% Root mean square of whisking
lenw = length(fwhisk);
wl = sr / 100;     % 10 ms windows (for 1 kHz samp. rate)
lle = floor(lenw/wl) * wl;
wh2 = reshape(fwhisk(1:lle),wl,lle/wl);
rms_whisk = sqrt(sum(wh2.^2)) / sqrt(wl);
rms_whisk = rms_whisk';
resp_phase2 = resp_phase(1:lle);
resp_phase2 = resp_phase2(round(wl/2):wl:end);   % downsampled resp. phase aligned to whisking RMS

% Respiration phase - whisking relationship
edges = -pi:2*pi/18:pi;     % phase histogram bin limits
cnts = (edges(1:end-1) + edges(2:end)) / 2;     % phase histogram bin centers
sfwhisk = fwhisk .^ 2;
spvr = [resp_phase sfwhisk];
spvr = sortrows(spvr,1);
sresp_phase = spvr(:,1);
ssfwhisk = spvr(:,2);    % phase-sorted whisking
[mn_phase mn_wh] = dhist(spvr,edges);
figure
plot([cnts cnts+2*pi],[mn_wh mn_wh],'r')      % conditional distribution: E(whisk|phi1<phase<phi2)
y_lim = ylim;
ylim([0 y_lim(2)])
[R p] = lincirc_corr2(mn_wh',cnts')

figure      % phase-sorted whisking
subplot(2,1,1)
plot(sresp_phase)
line(xlim,[0 0],'Color','g')
[b,bint,r,rint,stats] = regress(ssfwhisk,[ones(length(sresp_phase),1) ...
    sin(sresp_phase) cos(sresp_phase)]);      % correlation
R = sqrt(stats(1));
p = stats(3);
x_lim = xlim;
y_lim = ylim;
pl = oom(p);
text(x_lim(1)+(x_lim(2)-x_lim(1))*0.1,y_lim(1)+(y_lim(2)-y_lim(1))*0.95,...
    ['R = ' num2str(R) '  p ' pl],'Color','red')
subplot(2,1,2)
plot(ssfwhisk)

figure      % phase-sorted whisking RMS
spvr = [resp_phase2 rms_whisk];
spvr = sortrows(spvr,1);
sresp_phase = spvr(:,1);
ssfwhisk = spvr(:,2);    % phase-sorted whisking
subplot(2,1,1)
plot(sresp_phase)
line(xlim,[0 0],'Color','g')
[b,bint,r,rint,stats] = regress(ssfwhisk,[ones(length(sresp_phase),1) ...
    sin(sresp_phase) cos(sresp_phase)]);      % correlation
R = sqrt(stats(1));
p = stats(3);
x_lim = xlim;
y_lim = ylim;
pl = oom(p);
text(x_lim(1)+(x_lim(2)-x_lim(1))*0.1,y_lim(1)+(y_lim(2)-y_lim(1))*0.95,...
    ['R = ' num2str(R) '  p ' pl],'Color','red')
subplot(2,1,2)
plot(ssfwhisk)

figure
plot(resp_phase2,rms_whisk,'.')

figure      % smoothed version
[ssrp sds] = smooth(spvr(:,1),'circular',101);
subplot(2,1,1)
plot(ssrp)
hold on
plot(ssrp+sds,'Color',[204 204 204]/255)
plot(ssrp-sds,'Color',[204 204 204]/255)
line(xlim,[0 0],'Color','g')
subplot(2,1,2)
[ssfw sds] = smooth(spvr(:,2),'circular',101);
plot(ssfw)
hold on
plot(ssfw+sds,'Color',[204 204 204]/255)
plot(ssfw-sds,'Color',[204 204 204]/255)

% EEG phase
% eeg_phase = angle(hilbert(standardize(feeg)));    % Hilbert-transform
% figure
% plot(standardize(eeg))
% hold on
% plot(standardize(feeg),'r')
% plot(eeg_phase,'c')
% figure
% plot(standardize(whisk))
% hold on
% plot(eeg_phase,'r')

[eeg_phase inx] = lphase(eeg,sr,4,12,4);    % Hilbert-transform
figure
plot(standardize(eeg))
hold on
plot(standardize(feeg),'c')
plot(eeg_phase,'g')
plot(inx,eeg_phase(inx),'r')
figure
plot(standardize(whisk))
hold on
plot(eeg_phase,'r')

% EEG phase - whisking relationship
edges = -pi:2*pi/18:pi;     % phase histogram bin limits
cnts = (edges(1:end-1) + edges(2:end)) / 2;     % phase histogram bin centers
sfwhisk = fwhisk .^ 2;
eeg_phase2 = eeg_phase;
eeg_phase2(inx) = [];
sfwhisk2 = sfwhisk;
sfwhisk2(inx) = [];
inx2 = find(sfwhisk<10000);
eeg_phase2(inx2) = [];
sfwhisk2(inx2) = [];
spvr = [eeg_phase2 sfwhisk2];
spvr = sortrows(spvr,1);
sresp_phase = spvr(:,1);
ssfwhisk = spvr(:,2);    % phase-sorted whisking
[mn_phase mn_wh] = dhist(spvr,edges);
figure
plot([cnts cnts+2*pi],[mn_wh mn_wh],'r')      % conditional distribution: E(whisk|phi1<phase<phi2)
y_lim = ylim;
ylim([0 y_lim(2)])
[R p] = lincirc_corr2(mn_wh',cnts')

% -------------------------------------------------------------------------
function [mn_phase mn_rt] = dhist(spvr,edges)

n = length(edges);
mn_phase = zeros(1,n-1);
mn_rt = zeros(1,n-1);
for k = 2:n
    inx = find(spvr(:,1)>edges(k-1)&spvr(:,1)<edges(k));
    mn_phase(k-1) = length(inx);
    mn_rt(k-1) = mean(spvr(inx,2));
end

% -------------------------------------------------------------------------
function pl = oom(p)

if isequal(p,0)
    pl = '< 10^{-15}';
elseif p < 0.05 && p > 0.01
    pl = '< 0.05';
elseif p < 0.01 && p > 0.001
    pl = '< 0.01';
elseif p < 0.001 && p > 0.0001
    pl = '< 0.001';
elseif p < 0.0001
    ep = (-3:-1:-16);
    ords = 10 .^ ep;
    iop = find(ords>p,1,'last');
    epp = ep(iop);
    pl = ['< 10^{' num2str(epp) '}'];
else
    pl = ['= ' num2str(p)];
end

% -------------------------------------------------------------------------
function [ahee inx] = lphase(eeg,sr,fr1,fr2,mlt)
%LPHASE   EEG pahse.
%   [AHEE INX] = LPHASE(EEG,SR,FR1,FR2,MLT) returns time instances where
%   filtered EEG is below mean + MLT * SD threshold or the cycle length
%   falls out of the phase range (INX) along with EEG pjase series (AHEE).

% Filtering EEG
nqf = sr / 2;
flt = fir1(4096,[fr1 fr2]/nqf,'band');      % bandpass filtering
feeg = filtfilt(flt,1,eeg);
feeg = standardize(feeg);

% Hilbert transformation of the EEG
ahee = angle(hilbert(feeg));
aheedeg = ahee * (180 / pi);

% Check criteria:
% 1. discard cicles with EEG amp. lower then MLT * SD
% 2. discard cicles shorter then 1/fr2 s
% 3. discard cicles longer then 1/fr1 s
fn0 = valuecrossing(1:length(ahee),ahee',0,'down');
fn = round(fn0);
sd = std(feeg);
inx = [];
for k = 1:length(fn)-1
    seeg = feeg(fn(k)+1:fn(k+1));
    axs = max(seeg) - min(seeg);
    if (axs < mlt * sd)  || (fn(k+1) - fn(k) > 1 / fr1 * sr) || (fn(k+1) - fn(k) < 1 / fr2 * sr)
        inx = [inx fn(k)+1:fn(k+1)];
    end
end