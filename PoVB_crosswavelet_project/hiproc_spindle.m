%% DATAPATH

global DATAPATH

%% load Ivan 6

load([DATAPATH 'Hajni\EEGMPO_deltanorm\crosswavelet_spindle\Ivan06_11_5731_VPM_149_182_CWVARS.mat'])
cmx = MaxSpindlePower;
fname = 'Ivan6';
nucl = 'VPM';

%% load Ivan 9

load([DATAPATH 'Hajni\EEGMPO_deltanorm\crosswavelet_spindle\Ivan09_13_6163_VPM_195_256_CWVARS.mat'])
cmx1 = MaxSpindlePower;
load([DATAPATH 'Hajni\EEGMPO_deltanorm\crosswavelet_spindle\Ivan09_13_6163_VPM_274_364_CWVARS.mat'])
cmx2 = MaxSpindlePower;
cmx = [cmx1 cmx2];
fname = 'Ivan9';
nucl = 'VPM';

%% load Ivan 7

load([DATAPATH 'Hajni\EEGMPO_deltanorm\crosswavelet_spindle\Ivan07_09_5575_Po_366_501_CWVARS.mat'])
cmx = MaxSpindlePower;
fname = 'Ivan7';
nucl = 'Po';

%% load Ivan 32

load([DATAPATH 'Hajni\EEGMPO_deltanorm\crosswavelet_spindle\Ivan32_03_5215_Po_391_451_CWVARS.mat'])
cmx1 = MaxSpindlePower;
load([DATAPATH 'Hajni\EEGMPO_deltanorm\crosswavelet_spindle\Ivan32_03_5215_Po_554_610_CWVARS.mat'])
cmx2 = MaxSpindlePower;
cmx = [cmx1 cmx2];
fname = 'Ivan32';
nucl = 'Po';

%% load Ivan 38

load([DATAPATH 'Hajni\EEGMPO_deltanorm\crosswavelet_spindle\Ivan38_01_5289_Po_141_173_CWVARS.mat'])
cmx1 = MaxSpindlePower;
load([DATAPATH 'Hajni\EEGMPO_deltanorm\crosswavelet_spindle\Ivan38_01_5289_Po_657_697_CWVARS.mat'])
cmx2 = MaxSpindlePower;
load([DATAPATH 'Hajni\EEGMPO_deltanorm\crosswavelet_spindle\Ivan38_01_5289_Po_975_1080_CWVARS.mat'])
cmx3 = MaxSpindlePower;
load([DATAPATH 'Hajni\EEGMPO_deltanorm\crosswavelet_spindle\Ivan38_01_5289_Po_1143_1280_CWVARS.mat'])
cmx4 = MaxSpindlePower;
cmx = [cmx1 cmx2 cmx3 cmx4];
fname = 'Ivan38';
nucl = 'Po';

%% load Ivan 40

load([DATAPATH 'Hajni\EEGMPO_deltanorm\crosswavelet_spindle\Ivan40_01_5628_Po border cell_646_772_CWVARS.mat'])
cmx1 = MaxSpindlePower;
load([DATAPATH 'Hajni\EEGMPO_deltanorm\crosswavelet_spindle\Ivan40_01_5628_Po border cell_953_1121_CWVARS.mat'])
cmx2 = MaxSpindlePower;
load([DATAPATH 'Hajni\EEGMPO_deltanorm\crosswavelet_spindle\Ivan40_01_5628_Po border cell_1131_1346_CWVARS.mat'])
cmx3 = MaxSpindlePower;
load([DATAPATH 'Hajni\EEGMPO_deltanorm\crosswavelet_spindle\Ivan40_01_5628_Po border cell_1365_1480_CWVARS.mat'])
cmx4 = MaxSpindlePower;
cmx = [cmx1 cmx2 cmx3 cmx4];
fname = 'Ivan40';
nucl = 'PoVPM';

%% load Ivan 8

load([DATAPATH 'Hajni\EEGMPO_deltanorm\crosswavelet_spindle\Ivan08_cx_intra_1613_15_44_CWVARS.mat'])
cmx1 = MaxSpindlePower;
load([DATAPATH 'Hajni\EEGMPO_deltanorm\crosswavelet_spindle\Ivan08_cx_intra_1613_46_90_CWVARS.mat'])
cmx2 = MaxSpindlePower;
load([DATAPATH 'Hajni\EEGMPO_deltanorm\crosswavelet_spindle\Ivan08_cx_intra_1613_242_302_CWVARS.mat'])
cmx3 = MaxSpindlePower;
cmx = [cmx1 cmx2 cmx3];
fname = 'Ivan8';
nucl = 'Cx';

%% nov 23 c1

load([DATAPATH 'Hajni\EEGMPO_deltanorm\crosswavelet_spindle_newcells\23 nov 2005_c1__unit_S1_EEG_18_33s_CWVARS.mat'])
cmx1 = MaxSpindlePower;
load([DATAPATH 'Hajni\EEGMPO_deltanorm\crosswavelet_spindle_newcells\23 nov 2005_c1__unit_S1_EEG_39_56s_CWVARS.mat'])
cmx2 = MaxSpindlePower;
cmx = [cmx1 cmx2];
fname = 'nov 23 c1';
nucl = 'VPM';

%% nov 28

load([DATAPATH 'Hajni\EEGMPO_deltanorm\crosswavelet_spindle_newcells\28 nov 2005_c2__unit_S1_EEG_0_21s_CWVARS.mat'])
cmx = MaxSpindlePower;
fname = 'nov 28';
nucl = 'Po';

%% intra

load([DATAPATH 'Hajni\EEGMPO_deltanorm\crosswavelet_spindle_newcells\intra_4726__unit_S1_EEG_14_31s_CWVARS.mat'])
cmx1 = MaxSpindlePower;
load([DATAPATH 'Hajni\EEGMPO_deltanorm\crosswavelet_spindle_newcells\intra_4726__unit_S1_EEG_79_92s_CWVARS.mat'])
cmx2 = MaxSpindlePower;
cmx = [cmx1 cmx2];
fname = 'intra';
nucl = 'Po';

%% ivan 46

load([DATAPATH 'Hajni\EEGMPO_deltanorm\crosswavelet_spindle_newcells\Ivan_46_01_5736_unit__S1_EEG_845_881s_CWVARS.mat'])
cmx1 = MaxSpindlePower;
load([DATAPATH 'Hajni\EEGMPO_deltanorm\crosswavelet_spindle_newcells\Ivan_46_01_5736__unit_S1_EEG_993_1089_CWVARS.mat'])
cmx2 = MaxSpindlePower;
cmx = [cmx1 cmx2];
fname = 'ivan 46';
nucl = 'Po';

%% ivan 47

load([DATAPATH 'Hajni\EEGMPO_deltanorm\crosswavelet_spindle_newcells\Ivan_47_03_5762_unit__S1_EEG_124_200s_CWVARS.mat'])
cmx1 = MaxSpindlePower;
load([DATAPATH 'Hajni\EEGMPO_deltanorm\crosswavelet_spindle_newcells\Ivan_47_03_5762__unit_S1_EEG_222_244s_CWVARS.mat'])
cmx2 = MaxSpindlePower;
load([DATAPATH 'Hajni\EEGMPO_deltanorm\crosswavelet_spindle_newcells\Ivan_47_03_5762_unit__S1_EEG_374_424s_CWVARS.mat'])
cmx3 = MaxSpindlePower;
load([DATAPATH 'Hajni\EEGMPO_deltanorm\crosswavelet_spindle_newcells\Ivan_47_03_5762_unit__S1_EEG_547_606s_CWVARS.mat'])
cmx4 = MaxSpindlePower;
load([DATAPATH 'Hajni\EEGMPO_deltanorm\crosswavelet_spindle_newcells\Ivan_47_03_5762__unit_S1_EEG_617_693s_CWVARS.mat'])
cmx5 = MaxSpindlePower;
cmx = [cmx1 cmx2 cmx3 cmx4 cmx5];
fname = 'ivan 47';
nucl = 'Po';

%% trial

load([DATAPATH 'Hajni\EEGMPO_deltanorm\crosswavelet_spindle_newcells\trial 1_5697__unit_S1_EEG_0_54s_CWVARS.mat'])
cmx = MaxSpindlePower;
fname = 'trial';
nucl = 'Po';

%% xlswrite

global DATAPATH
xlsname = [DATAPATH 'Hajni\EEGMPO_deltanorm\crosswavelet_spindle_newcells\MaxSpindlePower.xls'];   % write results to excel file
if b_isfilename(xlsname)
    [ntz mtz atz] = xlsread(xlsname,'sheet1');
    pref = size(atz,1) + 1;
else
    pref = 1;
end
xlswrite(xlsname,{fname},'sheet1',['A' num2str(pref)])
xlswrite(xlsname,{nucl},'sheet1',['B' num2str(pref)])
xlswrite(xlsname,mean(cmx),'sheet1',['C' num2str(pref)])
xlswrite(xlsname,std(cmx),'sheet1',['D' num2str(pref)])
xlswrite(xlsname,std(cmx)/sqrt(length(cmx)),'sheet1',['E' num2str(pref)])