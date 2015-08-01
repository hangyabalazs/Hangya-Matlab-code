%%

eeg=CSC6_Samples(:);

%%

sr = 1/(CSC6_TimeStamps(2)-CSC6_TimeStamps(1))*1000000*512;
dt = 1 / sr;
time = (0:length(eeg)-1)*dt+CSC6_TimeStamps(1)/1000000;

%%

eeg = eeg(1:30:end);
sr = sr/30;
time = time(1:30:end);

%%

nqf = sr/2;      % filtering EEG
flt = fir1(4096,[20 35]/nqf,'band');      % lowpass filtering on 5 Hz
feeg = filtfilt(flt,1,eeg);

%%

nqf = sr/2;      % filtering EEG
flt = fir1(4096,[4 12]/nqf,'band');      % lowpass filtering on 5 Hz
feeg = filtfilt(flt,1,eeg);

%%

figure;
plot(time(2000000:end),eeg(2000000:end))
hold on
plot(time(2000000:end),feeg(2000000:end),'r')

%%

figure;
plot(time(1:1000000),eeg(1:1000000))
hold on
plot(time(1:1000000),feeg(1:1000000),'r')

%%

figure;
plot(time,eeg)
hold on
plot(time,feeg,'r')

%%

line([TS'; TS'],[zeros(length(TS),1)'; ones(length(TS),1)'*10000],'Color','c')

%%

aph = abs(hilbert(feeg));

figure;
plot(time,aph)
hold on
line([TS'; TS'],[zeros(length(TS),1)'; ones(length(TS),1)'*10000],'Color','c')

%%

[H1 H2 trsc] = somccg(TS,TS,1000);

%%

z = zeros(1,length(eeg));
z(round(TS)) = 1;
hold on
plot(time(1:1000000),z(1:1000000)*10000,'c')

%%

nqf = sr;      % filtering EEG
flt = fir1(4096,[0.5 4]/nqf,'band');      % lowpass filtering on 5 Hz
feeg = filtfilt(flt,1,eeg);

%%
beh_events=Events_TimeStamps(find(strcmp('Zoned Video: Zone1 Entered',Events_EventStrings)));

%%

vdisc = [];
for k = 1:length(beh_events)
    vdisc = [vdisc; TS(TS>beh_events(k)&TS<beh_events(k)+2)];
end

%%

[hang, hmvl, ftm, bahee] = somphase(eeg,TS,sr,time,3,8);

%%

[hang, hmvl, ftm, bahee] = somphase(eeg,TS,sr,time,20,35);

%%

x= sin(2*pi*28*[0:1/1000:50]).*(1-(sin(2*pi*0.5*[0:1/200:50])/10));
x= sin(2*pi*28*[0:1/1000:50]);
figure;plot(x)

nqf = 200/2;      % filtering EEG
flt = fir1(512,[20 35]/nqf,'band');      % lowpass filtering on 5 Hz
fx = filtfilt(flt,1,x);

time = 0:1/1000:50;

figure;
plot(time,x)
hold on
plot(time,fx,'r')
vd = 1/540:1/540:50;
z = zeros(1,length(x));
z(ceil(vd)) = 1;
hold on
plot(time,z,'c')


[hang, hmvl, ftm, bahee] = somphase(x',vd',1000,time,20,35);

%%

figure
[mu,kappa,Value,p,Rsquare] = b_watson2(bahee,1)