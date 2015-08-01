%% import

unitr = data(:,2);

figure
plot(unitr)

%% filter

sr = 10000;
nqf = sr / 2;
fl = fir1(512,200/sr,'high');
unit = filtfilt(fl,1,unitr)';
fl = fir1(512,200/sr,'low');
eeg = filtfilt(fl,1,unitr)';

figure
plot(eeg)
figure
plot(unit)

%% discrimination

vdisc = disc(unit,0.7);

%% ISI histogram

isi = diff(vdisc);
figure
hist(isi(isi<1*sr),1000)
figure
hist(isi,1000)

%% instant. freq.

instfrek = [];lenu=length(eeg);
isi = diff(vdisc);
for i = 1:length(vdisc)-1
    instfrek(vdisc(i):vdisc(i+1)) = 1 / isi(i);
end
instfrek(1:vdisc(1)-1) = 1 / vdisc(1);
instfrek(vdisc(end):lenu) = 1 / (lenu - vdisc(end));
instfrek2 = instfrek * 10000;
figure;plot(instfrek2)

%% wavelet

lenu = length(eeg);
[powunit,phaseunit,f] = unitwavelet(vdisc,lenu,sr);
figure
imagesc(powunit)
b_rescaleaxis('Y',f)
[poweeg,phaseeeg,f] = eegwavelet(eeg(1:10:end),sr/10);
figure
imagesc(poweeg)
b_rescaleaxis('Y',f)
figure
imagesc(poweeg(:,62*sr/10:90*sr/10))
b_rescaleaxis('Y',f)
figure
imagesc(poweeg,[0 230])
b_rescaleaxis('Y',f)

%% unit theta/delta power

pwind1 = find(f<6,1,'first');        % theta band
pwind2 = find(f>2.5,1,'last');
pwind3 = find(f<2.5,1,'first');      % delta band
pwind4 = find(f>0.5,1,'last');
tpd = sum(powunit(pwind1:pwind2,:)) ./ sum(powunit(pwind3:pwind4,:));
figure
plot(tpd)

%% EEG theta/delta power

pwind1 = find(f<6,1,'first');        % theta band
pwind2 = find(f>2.5,1,'last');
pwind3 = find(f<2.5,1,'first');      % delta band
pwind4 = find(f>0.5,1,'last');
tpd = sum(poweeg(pwind1:pwind2,:)) ./ sum(poweeg(pwind3:pwind4,:));
figure
plot(tpd)

%% phase

eeg_theta = eeg(60*sr:90*sr);
vdisc_theta = vdisc(vdisc>60*sr&vdisc<90*sr)-60*sr;

% Filtering EEG
flt = fir1(512,[2/nqf 8/nqf]);
feeg = filtfilt(flt,1,eeg_theta);

% Hilbert transformation of the eeg
ahee = angle(hilbert(feeg));
aheedeg = ahee * (180 / pi);

% Phase angles - Hilbert
bahee = ahee(vdisc_theta(vdisc_theta>0&vdisc_theta<length(eeg_theta)));
n = length(bahee);

ftm = sum(exp(1).^(i*bahee)) / n;    % first trigonometric moment
hang = angle(ftm);   % mean angle
hmvl = abs(ftm);     % mean resultant length

figure
[nm xout] = hist(bahee,20);
bar([xout xout+2*pi],[nm nm])

%% firing rate

seglen = 5 * sr;        % 5 sec. long segments
len = length(unit);
olp = 5;    % 80% overlapping windows
ind1 = [1:seglen/olp:len-seglen+1];
ind2 = ind1 + seglen - 1;
lenr = length(ind1);
frate = zeros(1,lenr);
for k = 1:lenr
    vd = vdisc(vdisc>ind1(k)&vdisc<ind2(k)) - ind1(k);
    loceeg = eeg(ind1(k):ind2(k));
    efflen = (vd(end) - vd(1)) / sr;
    frate(k) = length(vd) / efflen;
end
figure
plot(frate)
ylabel('firing rate')

%% sharp waves

flt = fir1(2048,[90/nqf 140/nqf]);
feeg = filtfilt(flt,1,eeg);
for k = 1:size(SharpWaveSegments,2)
    figure
    subplot(2,1,1)
    loceeg = feeg(SharpWaveSegments(1,k)-1000:SharpWaveSegments(2,k)+1000);
    plot(loceeg)
    subplot(2,1,2)
    locunit = unit(SharpWaveSegments(1,k)-1000:SharpWaveSegments(2,k)+1000);
    plot(locunit)
end
figure
wn = 10000;
z = zeros(1,2*wn);
for k = 1:size(SharpWaveSegments,2)
    cnt = round(mean([SharpWaveSegments(1,k),SharpWaveSegments(2,k)]))
    loceeg = feeg(cnt-wn:cnt+wn);
    subplot(3,1,1)
    hold on
    plot(loceeg)
    subplot(3,1,2)
    hold on
    locunit = unit(cnt-wn:cnt+wn);
    plot(locunit)
    locvdisc = (vdisc(vdisc>cnt-wn&vdisc<cnt+wn)) - cnt + wn;
    z(locvdisc) = 1;
end
z2 = reshape(z,2*wn/20,20);
subplot(3,1,3)
bar(1:20,sum(z2))

figure
wn = 10000;
z = zeros(1,2*wn);
for k = 1:size(SharpWaveSegments,2)
    cnt = SharpWaveSegments(1,k)
    loceeg = feeg(cnt-wn:cnt+wn);
    subplot(3,1,1)
    hold on
    plot(loceeg)
    subplot(3,1,2)
    hold on
    locunit = unit(cnt-wn:cnt+wn);
    plot(locunit)
    locvdisc = (vdisc(vdisc>cnt-wn&vdisc<cnt+wn)) - cnt + wn;
    z(locvdisc) = 1;
end
z2 = reshape(z,2*wn/20,20);
subplot(3,1,3)
bar(1:20,sum(z2))

figure
wn = 10000;
z = zeros(1,2*wn);
for k = 1:size(SharpWaveSegments,2)
    cnt = SharpWaveSegments(2,k);
    loceeg = feeg(cnt-wn:cnt+wn);
    subplot(3,1,1)
    hold on
    plot(loceeg)
    subplot(3,1,2)
    hold on
    locunit = unit(cnt-wn:cnt+wn);
    plot(locunit)
    locvdisc = (vdisc(vdisc>cnt-wn&vdisc<cnt+wn)) - cnt + wn;
    z(locvdisc) = 1;
end
z2 = reshape(z,2*wn/20,20);
subplot(3,1,3)
bar(1:20,sum(z2))
