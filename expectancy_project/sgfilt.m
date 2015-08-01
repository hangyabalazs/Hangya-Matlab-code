%% filtertest

x2=sin((1:2000)/2000*2*pi*4);
X=[zeros(1,4000) x2 x2]+randn(1,8000)/10;
figure;plot(X)
flt=fir1(512,[0.5 12]/1000);
fX=filtfilt(flt,1,X);
figure;plot(fX)
X=[zeros(1,4000) x2 x2]+randn(1,8000)/10;
fX=filtfilt(flt,1,X);
figure;plot(fX)
axis([3000 4000 0.5 -0.5])
axis([3000 4000 -0.1 0.1])

%% test2

x2 = sin((1:2000)/2000*2*pi*4);     % sampling rate: 2000 Hz
flt = fir1(512,[0.5 12]/1000);
fX = zeros(100,8000);
for k = 1:100
    X = [zeros(1,4000) x2 x2]+randn(1,8000)/10;     % 8 sec.
    fX(k,:) = filtfilt(flt,1,X);
end
fXX = mean(fX);
figure
plot(fXX)

%% test3

x2 = sin((1:2000)/2000*2*pi*4) * 50;     % sampling rate: 2000 Hz
flt = fir1(512,[0.5 12]/1000);
fX = zeros(100,8000);
for k = 1:100
    X = [zeros(1,4000) x2 x2]+randn(1,8000)/10;     % 8 sec.
    fX(k,:) = filtfilt(flt,1,X);
end
fXX = mean(fX);
figure
plot(fXX)

%% test4

x2 = sin((1:2000)/2000*2*pi*4) * 50;     % sampling rate: 2000 Hz
flt = fir1(512,[0.5 12]/1000);
fX = zeros(100,8000);
for k = 1:100
    X = [zeros(1,4000) x2 x2];     % 8 sec.
    fX(k,:) = filtfilt(flt,1,X);
end
fXX = mean(fX);
figure
plot(fXX)

%% test5

n = normpdf((1:4000),2000,1000);
x2 = sin((1:2000)/2000*2*pi*4) * 50;     % sampling rate: 2000 Hz
flt = fir1(512,[0.5 12]/1000);
fX = zeros(100,8000);
for k = 1:100
    X = [zeros(1,4000) [x2 x2].*n]+randn(1,8000)/10;     % 8 sec.
    fX(k,:) = filtfilt(flt,1,X);
end
fXX = mean(fX);
figure
plot(fXX)