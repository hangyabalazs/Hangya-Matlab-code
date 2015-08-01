function snwphase3

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
% flt = fir1(1024,[4 12]/nqf,'band');     % filter EEG in the theta band
% feeg = filtfilt(flt,1,eeg);
% figure
% plot(whisk)
% hold on
% plot(fwhisk,'r')
% figure
% plot(eeg)
% hold on
% plot(feeg,'r')

% Respiration phase
resp_phase = angle(hilbert(standardize(fsniff)));    % Hilbert-transform
% figure
% plot((sniff-mean(sniff))/std(sniff))
% hold on
% plot((fsniff-mean(fsniff))/std(fsniff),'r')
% plot(resp_phase,'c')
% figure
% plot((sniff-mean(sniff))/std(sniff))
% hold on
% plot(resp_phase,'r')

% % Discriminate whisking
% vdisc = b_udisc(fwhisk');
% 
% % 'Burst detection'
% isi = diff(vdisc);
% figure;
% lni = 900;
% plot(isi(1:lni-1),isi(2:lni),'.')

% Root mean square of whisking
lenw = length(fwhisk);
wl = sr / 100;     % 10 ms windows (for 1 kHz samp. rate)
lle = floor(lenw/wl) * wl;
wh2 = reshape(fwhisk(1:lle),wl,lle/wl);
rms_whisk = sqrt(sum(wh2.^2)) / sqrt(wl);
rms_whisk = rms_whisk';
% flt = fir1(512,10/(nqf/wl),'high');
% rms_whisk = filtfilt(flt,1,rms_whisk);  % higi-pass filter at 10 Hz
resp_phase2 = resp_phase(1:lle);
resp_phase2 = resp_phase2(round(wl/2):wl:end);   % downsampled resp. phase aligned to whisking RMS

% Detect whisking
[st nd pk fwhisk_] = sspw3(whisk,sr,10,10,inf,4,1);
figure
plot(rms_whisk)
hold on
line([pk/wl; pk/wl],[zeros(size(pk)); ones(size(pk))*10000],'Color','green')

% Get resp. cycles
fn0 = valuecrossing(1:length(resp_phase2),resp_phase2',0,'down');
fn = round(fn0);
figure
line([fn; fn],[ones(size(fn))-4; ones(size(fn))+3],'Color','green')
hold on
plot(resp_phase2)

% Sort whisking events
vdisc = round(pk/wl);
figure
line([vdisc; vdisc],[zeros(size(pk)); ones(size(pk))*2000],'Color','green')
hold on
plot(fsniff(1:wl:end))
grp = zeros(size(vdisc));
for k = 1:length(vdisc)
    ffv = find(fn<vdisc(k),1,'last');
    fn_4pre = fn(ffv-3);
    fn_pppre = fn(ffv-2);
    fn_ppre = fn(ffv-1);
    fn_pre = fn(ffv);
    fn_post = fn(ffv+1);
    fn_ppost = fn(ffv+2);
    fn_pppost = fn(ffv+3);
    fn_4post = fn(ffv+4);
    lv_pppre = length(vdisc(vdisc>fn_4pre&vdisc<fn_pppre));
    lv_ppre = length(vdisc(vdisc>fn_pppre&vdisc<fn_ppre));
    lv_pre = length(vdisc(vdisc>fn_ppre&vdisc<fn_pre));
    lv = length(vdisc(vdisc>fn_pre&vdisc<fn_post));
    lv_post = length(vdisc(vdisc>fn_post&vdisc<fn_ppost));
    lv_ppost = length(vdisc(vdisc>fn_ppost&vdisc<fn_pppost));
    lv_pppost = length(vdisc(vdisc>fn_pppost&vdisc<fn_4post));
    switch lv
        case 2
            grp(k) = 2;
        case 3
            grp(k) = 3;
        case 1
            if lv_pre >= 1 || lv_post >= 1
                grp(k) = (1);
            else
                if lv_ppre >= 1 || lv_ppost >= 1
                    grp(k) = 0.5;
                else
                    if lv_pppre >= 1 || lv_pppost >= 1
                        grp(k) = 1/3;
                    end
                end
            end
    end
end

