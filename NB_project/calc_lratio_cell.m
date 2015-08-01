%%

k=1
[DM LRC VALID_CHANNELS] = LRatio2(CELLIDLIST{k},'feature_names',{'WavePC1' 'Energy'},'valid_channels',[0 1 1 1])
setvalue(CELLIDLIST{k},'ID_PC',DM)
setvalue(CELLIDLIST{k},'Lr_PC',LRC)
[DM LRC VALID_CHANNELS] = LRatio2(CELLIDLIST{k},'feature_names',{'Amplitude' 'Energy'},'valid_channels',[0 1 1 1])
setvalue(CELLIDLIST{k},'Lr_amp',LRC)
setvalue(CELLIDLIST{k},'ID_amp',DM)