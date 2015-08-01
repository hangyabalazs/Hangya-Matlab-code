%% spike width distribution

ss = getvalue('SpikeShape');
ss2 = nancell2struct(ss);
ss2 = nancell2struct2(ss);
sw = ss2.PeakToPostValleyTime;
sw = [ss2.PeakToPostValleyTime];
figure;hist(sw,22)

%% cells in behavior

selstr = ['"ChAT+"==1&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("session_type",{''gonogo'',''feedbackdelay''})&' ...
    'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})'];
ChAT = selectcell(selstr);   % cell IDs for ChAT cells
ChAT = [ChAT 'n045_121217x_4.6'];   % clustered based on light-evoked spikes
selstr = ['"pChAT+"==1&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("session_type",{''gonogo'',''feedbackdelay''})&' ...
    'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})'];
pChAT = selectcell(selstr);  % putative
allChAT = [ChAT pChAT];
selstr = ['"ChAT+"==0&"pChAT+"==0&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("session_type",{''gonogo'',''feedbackdelay''})&' ...
    'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})'];
NT = selectcell(selstr);   % cell IDs for non-tagged cells

%% spike width per group

ss_allChAT = getvalue('SpikeShape',[ChAT pChAT]);
ss2_allChAT = nancell2struct2(ss);
SpikeWidth_allChAT = ss2.PeakToPostValleyTime;
sw = [ss2.PeakToPostValleyTime];

%% spike width per group

SpikeShape_NT = getvalue('SpikeShape',NT);
SpikeShape_NT = nancell2struct(SpikeShape_NT);
SpikeWidth_NT = [SpikeShape_NT.PeakToPostValleyTime];
Spike_NT = {SpikeShape_NT.Spike};
Spike_NT = cell2mat(cellfun(@(s)zscore(s'),Spike_NT,'UniformOutput',false))';

SpikeShape_allChAT = getvalue('SpikeShape',allChAT);
SpikeShape_allChAT = nancell2struct(SpikeShape_allChAT);
SpikeWidth_allChAT = [SpikeShape_allChAT.PeakToPostValleyTime];
Spike_allChAT = {SpikeShape_allChAT.Spike};
Spike_allChAT = cell2mat(cellfun(@(s)zscore(s'),Spike_allChAT,'UniformOutput',false))';

SpikeShape_ChAT = getvalue('SpikeShape',ChAT);
SpikeShape_ChAT = nancell2struct(SpikeShape_ChAT);
SpikeWidth_ChAT = [SpikeShape_ChAT.PeakToPostValleyTime];
Spike_ChAT = {SpikeShape_ChAT.Spike};
Spike_ChAT = cell2mat(cellfun(@(s)zscore(s'),Spike_ChAT,'UniformOutput',false))';

SpikeShape_pChAT = getvalue('SpikeShape',pChAT);
SpikeShape_pChAT = nancell2struct(SpikeShape_pChAT);
SpikeWidth_pChAT = [SpikeShape_pChAT.PeakToPostValleyTime];
Spike_pChAT = {SpikeShape_pChAT.Spike};
Spike_pChAT = cell2mat(cellfun(@(s)zscore(s'),Spike_pChAT,'UniformOutput',false))';

% SpikeShape_PV = getvalue('SpikeShape',PV);
% SpikeShape_PV = nancell2struct(SpikeShape_PV);
% SpikeWidth_PV = [SpikeShape_PV.PeakToPostValleyTime];
% Spike_PV = {SpikeShape_PV.Spike};
% Spike_PV = cell2mat(cellfun(@(s)zscore(s'),Spike_PV,'UniformOutput',false))';

SpikeWidth_allcells = [SpikeWidth_NT SpikeWidth_allChAT];

%% FR

Baseline_NT = getvalue('Baseline',NT);
Baseline_ChAT = getvalue('Baseline',ChAT);
Baseline_pChAT = getvalue('Baseline',pChAT);
Baseline_allChAT = [Baseline_ChAT; Baseline_pChAT];
Baseline_allcells = [Baseline_NT; Baseline_allChAT];

%% target population: fast NS cells

NS = NT(Baseline_NT>15&SpikeWidth_NT'<200);


%% tuning curves

nbtuningcurves3(NS,true)