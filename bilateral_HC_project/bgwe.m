function bgwe(data,trials)
%BGWE   Wavelet entropy on Buzsaki's EEG files.
%   BGWE calculates wavelet, FFT and filtered signals for low and high
%   gamma band.
%
%   See also BGWE2.

% Get trial EEG
lt = length(trials);
lt = 1;
for k = lt
    tS = round(trials{k}.tStart*4/5);  % indices are transformed to match new samp. rate
    tE = round(trials{k}.tEnd*4/5);
    eeg = data(4,tS:tE);
    
% Filter
    sr = 1000;
    nqf = sr / 2;
    flt = fir1(1024,[30 40]/nqf,'band');      % bandpass filtering
    feeg_low = filtfilt(flt,1,eeg);
    feeg_low = (feeg_low - mean(feeg_low)) / std(feeg_low);
    flt = fir1(1024,[80 100]/nqf,'band');      % bandpass filtering
    feeg_high = filtfilt(flt,1,eeg);
    feeg_high = (feeg_high - mean(feeg_high)) / std(feeg_high);
    
% Wavelet
    [wavea_abs wavea_phase f] = eeg_wavelet(eeg); % eeg sampled on 1000 Hz
    figure
    imagesc(wavea_abs)
    b_rescaleaxis('Y',f)
    
    fnd = find(f>30);    % frequency band bounderies
    pwind1 = fnd(end);
    fnd = find(f<100);
    pwind2 = fnd(1);
    figure
    S1 = subplot(3,1,1);
    imagesc(wavea_abs(1:pwind1,:))
    b_rescaleaxis('Y',f)
    S2 = subplot(3,1,2);
    plot(eeg)
    hold on
    plot(feeg_low*1000,'Color',[204 204 204]/255)
    x_lim = xlim;
    xlim([x_lim(1) x_lim(1)+length(eeg)])
    S2 = subplot(3,1,3);
    plot(eeg)
    hold on
    plot(feeg_high*1000,'Color',[204 204 204]/255)
    x_lim = xlim;
    xlim([x_lim(1) x_lim(1)+length(eeg)])
    
% FFT
    eeg2 = eeg(1:5:end);        % downsample on 200 Hz
    [y,w] = b_fft2(eeg2,200);     % FFT
    inxs_low = w > 30 & w < 40;
    inxs_high = w > 80 & w < 100;
    inxs2 = w > 5 & w < 100;
    integrated_power_low = sum(y(inxs_low));
    integrated_power_high = sum(y(inxs_high));
    H = figure;
    plot(w(inxs2),y(inxs2))
    x_lim = xlim;
    y_lim = ylim;
    text((x_lim(2)-x_lim(1))/2+x_lim(1),(y_lim(2)-y_lim(1))*0.75+y_lim(1),...
        ['30-40 Hz: ' num2str(integrated_power_low/1000000) ' * 10^6'])
    text((x_lim(2)-x_lim(1))/2+x_lim(1),(y_lim(2)-y_lim(1))*0.70+y_lim(1),...
        ['80-100 Hz: ' num2str(integrated_power_high/1000000) ' * 10^6'])
end



% -------------------------------------------------------------------------
% WAVELET
% -------------------------------------------------------------------------
function [pow,phase,f] = eeg_wavelet(dat)

% Prepare for wavelet transformation
variance = std(dat) ^ 2;
dat = (dat - mean(dat)) / sqrt(variance) ;
n = length(dat);
dt = 1 / 1000;
pad = 1;
dj = 0.08;    
j1 = ceil((1/dj) * log2(n/2));
j1 = ceil(j1);
j = (0:j1);
s0 = 2 * dt; 
s = s0 .* 2 .^ (j * dj);
omega0 = 6;
c = 4 * pi / (omega0 + sqrt(2+omega0^2));
fperiod = c .* s;
f = 1 ./ fperiod;
lag1 = 0.72;
param = -1;
mif = 0.5;          %minimal intresting frequency
mis = find(f>mif);
mis = mis(end);     %maximal intristing scale
mother = 'Morlet';

% Wavelet transformation
[wave,period,scale,coi] = b_wavelet_new3(dat,dt,pad,dj,s0,j1,mother,param,mis);
pow = abs(wave).^2;
phase = angle(wave);