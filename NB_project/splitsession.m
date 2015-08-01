%%

ct = 4*10^7;
inx = find(FeatureTimestamps>ct);
FeatureTimestamps = FeatureTimestamps(inx);
FeatureData = FeatureData(inx,:);
FeatureIndex = FeatureIndex(inx);

%% n026

save('c:\Balazs\_analysis\NB\junk\120222a\TT6_Amplitude.mat')

%%

save('c:\Balazs\_analysis\NB\junk\120222a\TT6_Energy.mat')

%%

save('c:\Balazs\_analysis\NB\junk\120222a\TT6_Time.mat')

%%

save('c:\Balazs\_analysis\NB\junk\120222a\TT6_wavePC1.mat')

%% step 1

fn = ['c:\Balazs\_analysis\NB\junk\120222a - Copy\TT6.ntt'];
[TimeStamp, ChanNum, SampleFrequency, NumValSamples, Samples, NlxHeader] = ...
    Nlx2MatSpike(fn,[1 1 1 1 1],1,1,1);

%% step 2

ct = 5*10^9;    % MClust timestamp * 100
inx = find(TimeStamp>ct);
Samples = Samples(:,:,inx);
TimeStamp = TimeStamp(inx);
SampleFrequency = SampleFrequency(inx);
ChanNum = ChanNum(inx);
NumValSamples = NumValSamples(:,inx);

%% step 3

fn = 'c:\Balazs\_analysis\NB\junk\120222a - Copy\TT6_.ntt';
Mat2NlxTT(fn,0,1,1,length(TimeStamp),[1 1 1 1 1 1],...
    TimeStamp,ChanNum,SampleFrequency,NumValSamples,Samples,NlxHeader)

%%

fn = ['f:\NB_cellbase\n029\120204b - Copy\TT3_.ntt'];
[TimeStamp, ChanNum, SampleFrequency, NumValSamples, Samples, NlxHeader] = ...
    Nlx2MatSpike(fn,[1 1 1 1 1],1,1,1);