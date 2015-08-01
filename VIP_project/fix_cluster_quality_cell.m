%%

cellid = 'ml01_130509b_2.10';
[ID_PC Lr_PC] = Lratio(cellid,{'Amplitude' 'WavePC1'})
[ID_amp Lr_amp] = Lratio(cellid,{'Amplitude' 'Energy'})
setvalue(cellid,'ID_PC',ID_PC)
setvalue(cellid,'Lr_PC',Lr_PC)
setvalue(cellid,'Lr_amp',Lr_amp)
setvalue(cellid,'ID_amp',ID_amp)
num2cell(listcell(cellid))