%% load

cd('c:\Balazs\_data\Intan\nba47\130527a')
[t,amps,data,aux] = read_intan_data;
cd('C:\MATLAB_R2010a\work\Balazs\Intan')
dt = data(:,4);

%% filter

% LFP
sr = 25000;
nqf = sr / 2;
flt = fir1(512,[0.1 600]/nqf,'band');
lfp = filtfilt(flt,1,dt);
lfp2 = lfp(1:25:end);   % downsample on 1000 Hz

% Unit
flt = fir1(2^7,[400 7000]/nqf,'band');
unit = filter(flt,1,dt);

%%

[b,a] = butter(3,[400 7000]/nqf,'bandpass');   % Butterworth filter
unit1 = filter(b,a,dt);


%% plot

% Plot LFP
figure;
plot(lfp2(1:100000/25))

% Plot unit
figure;
plot(unit(1:100000));