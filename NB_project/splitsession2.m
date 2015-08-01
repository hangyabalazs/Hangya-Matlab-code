%% step 1

fn = ['c:\Balazs\_analysis\NB\junk\121229g - Copy\TT1.ntt']; % n046
[TimeStamp, ChanNum, SampleFrequency, NumValSamples, Samples, NlxHeader] = ...
    Nlx2MatSpike(fn,[1 1 1 1 1],1,1,1);

%% step 2

ct = 3.4*10^10;    % MClust timestamp * 100
inx = find(TimeStamp<ct);
Samples = Samples(:,:,inx);
TimeStamp = TimeStamp(inx);
SampleFrequency = SampleFrequency(inx);
ChanNum = ChanNum(inx);
NumValSamples = NumValSamples(:,inx);

%% step 3

fn = 'c:\Balazs\_analysis\NB\junk\121229g - Copy\TT1_.ntt';
Mat2NlxTT(fn,0,1,1,length(TimeStamp),[1 1 1 1 1 1],...
    TimeStamp,ChanNum,SampleFrequency,NumValSamples,Samples,NlxHeader)

%% TrialEvents

TE = load('C:\Balazs\_analysis\NB\junk\121229g - Copy\TrialEvents.mat');
ctt = ct / 1000000;
tinx = find(TE.TrialStart<ctt);
fld = fieldnames(TE);
for k = 1:length(fld)
    TE.(fld{k}) = TE.(fld{k})(tinx);
end
save('C:\Balazs\_analysis\NB\junk\121229g - Copy\TrialEvents.mat','-struct','TE')