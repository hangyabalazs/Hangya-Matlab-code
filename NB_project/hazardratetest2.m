%% load

TE = load('C:\Balazs\_data\NB\NB_callbase\n013\110706a\TE.mat');
ITIs = TE.TotalITI;
figure;plot(TE.TotalITI(TE.TotalITI<10),TE.GoRT(TE.TotalITI<10),'.')

%% Failure distribution

ITIMin = min(ITIs);
ITIMax = max(ITIs);
NumTrials = length(ITIs);

%% Failure density function

dt = 0.5;
times = ITIMin-dt:dt:ITIMax+dt;
cnts = (times(1:end-1) + times(2:end)) / 2;
ft = histc(ITIs,times);
ft = ft(1:end-1);
ft = ft / sum(ft);
figure
bar(cnts,ft)

%% Hazard rate

Ft = cumsum(ft);
Rt = 1 - Ft;
ht = (Rt(1:end-1) - Rt(2:end)) ./ (dt * Rt(1:end-1));
figure
plot(cnts(1:end-1),ht)