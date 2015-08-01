%% load

TE = load('C:\Balazs\_data\NB\NB_callbase\n013\110628a\TrialEvents.mat');
load('C:\Balazs\_data\NB\NB_callbase\n013\110628a\EVENTSPIKES5_4.mat')
spikes_stimon = event_stimes{1};

%% prestimulus frequency

NUMtrials = length(spikes_stimon);
prestimfreq = nan(1,NUMtrials);
for k = 1:NUMtrials
    lspikes = spikes_stimon{k};
    lspikes2 = lspikes(lspikes>-2&lspikes<0);   % time window: one sec before stimulus onset
    prestimfreq(k) = length(lspikes2) / 2;
end

ishit = logical(nan2zero(TE.Hit));
ismiss = logical(nan2zero(TE.Miss));
isfa = logical(nan2zero(TE.FalseAlarm));
iscr = logical(nan2zero(TE.CorrectRejection));

%% conditional distributions

figure
[nm xout] = hist(prestimfreq(isfa));
plot(xout,nm)
hold on
[nm xout] = hist(prestimfreq(iscr));
plot(xout,nm,'r')

%% boxplot

H = figure;
b1 = prestimfreq(iscr);
b2 = prestimfreq(isfa);
boxplot([b1 b2],[zeros(size(b1)) ones(size(b2))],'labels',[{'CR'} {'FA'}]);
b_ranksum2(prestimfreq(iscr),prestimfreq(isfa))

%% firing rate -  go RT

is60 = TE.StimulusDuration == 30;

RT = TE.GoRT;
HitRT = RT(ishit&is60);
HitFR = prestimfreq(ishit&is60);
figure;
plot(HitFR,HitRT,'.')

%% firing rate -  no-go RT

is60 = TE.StimulusDuration == 30;

NoGoRT = TE.NoGoRT;
FART = NoGoRT(isfa&is60);
FAFR = prestimfreq(isfa&is60);
figure;
plot(FAFR,FART,'.')