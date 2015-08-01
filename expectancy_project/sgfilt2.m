%% test filter on orig. data

dim1 = 4;       % 91%
dim2 = 1;       % Fz
dat = squeeze(data(dim1,dim2,:,:))';
sr = 1000;
nqf = sr / 2;
fraw = zeros(100,3000);
nraw = zeros(100,3000);
flt = fir1(512,[0.5 12]/nqf);
for k = 1:100
    raw = dat(k,:);
    m(k) = raw(1500);
    raw(1:1500) = randn(1,1500)+raw(1500);
    nraw(k,:) = raw;
    fraw(k,:) = filtfilt(flt,1,raw);
end
frr = mean(fraw);
figure
plot(mean(nraw))
hold on
plot(frr,'r')

%% plot mean across trials

mraw = mean(dat);
sr = 1000;
nqf = sr / 2;
flt = fir1(512,[0.5 12]/nqf);
fmr = filtfilt(flt,1,mraw);

figure
plot(mraw)
hold on
plot(fmr,'r')

%% plot mean filterd data

mraw = mean(dat);
sr = 1000;
nqf = sr / 2;
flt = fir1(512,[0.5 12]/nqf);
fmr = zeros(100,3000);
for k = 1:100
    fmr(k,:) = filtfilt(flt,1,dat(k,:));
end
figure
plot(mraw)
hold on
plot(mean(fmr),'r')

%% test with SG's filter

dim1 = 4;       % 91%
dim2 = 1;       % Fz
dat = squeeze(data(dim1,dim2,:,:))';
sr = 1000;
Fs = sr;
Fp1 = 4;
Fp2 = 7;
nqf = sr / 2;
fraw = zeros(100,3000);
nraw = zeros(100,3000);
flt = fir1(512,[0.5 12]/nqf);
for k = 1:100
    raw = dat(k,:);
    m(k) = raw(1500);
%     raw(1:1500) = randn(1,1500) + raw(1500);
    nraw(k,:) = raw;
    fraw(k,:) = bandpassFilter(raw,1000,4,7);
end
frr = mean(fraw);
figure
plot(mean(nraw))
hold on
plot(frr,'r')

%% test with fir1 filter (4-7 Hz)

dim1 = 4;       % 91%
dim2 = 1;       % Fz
dat = squeeze(data(dim1,dim2,:,:))';
sr = 1000;
Fs = sr;
Fp1 = 4;
Fp2 = 7;
nqf = sr / 2;
flt = fir1(512,[4 7]/nqf);
fraw = zeros(100,3000);
nraw = zeros(100,3000);
flt = fir1(512,[0.5 12]/nqf);
for k = 1:100
    raw = dat(k,:);
    m(k) = raw(1500);
%     raw(1:1500) = randn(1,1500) + raw(1500);
    nraw(k,:) = raw;
    fraw(k,:) = filtfilt(flt,1,raw);
end
frr = mean(fraw);
figure
plot(mean(nraw))
hold on
plot(frr,'r')