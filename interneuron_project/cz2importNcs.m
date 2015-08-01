%%

[Timestamp, ChanNum, SampleFrequency, NumValSamples, Samples, NlxHeader] = ...
    Nlx2MatCSC('X:\In_Vivo\balazs\_analysis\Czurko2\INTERNEURONS\_osc_int_EEG\acin11s031g\CSC09.Ncs',[1 1 1 1 1],1,1,1);
Samples2=Samples(:);