%% to get reliability stats

choosecb('NB')
selstr = ['"ChAT+"==1&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("session_type",{''gonogo'',''feedbackdelay''})&' ...
    'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})'];
ChAT = selectcell(selstr);   % cell IDs for ChAT cells (n = 22)
ChAT = [ChAT 'n045_121217x_4.6'];   % clustered based on light-evoked spikes

global DATAPATH
impdir = fullfile([DATAPATH,'NB\taggingsummary_newdata\reliability2\mat_allcells\']);
dr = dir(impdir);
files = {dr(3:end).name};
NumCells = length(ChAT);
stim_freq = cell(1,NumCells);
Rs = nan(NumCells,4);
for iC = 1:NumCells
    cellid = ChAT{iC};
    if isequal(cellid,'n071_141218a_6.1')   % choose session w multiple stim frequencies
        cellid = 'n071_141229a_6.2';
    end
    cellidt = regexprep(cellid,'\.','_');
    fnm = fullfile(impdir,[cellidt '_tagging_reliability.mat']);
    load(fnm)
%     figure
%     plot(burst_types/2,Reliability)
    disp([num2str(iC) '  ' cellid])
%     keyboard
    stim_freq{iC} = burst_types / 2;
%     burst_types / 2
    Rs(iC,1:length(Reliability)) = Reliability;
end

choosecb('HDB')
ChAT = selectcell(['"ChAT+"==1&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("session_type",{''gonogo'',''feedbackdelay''})&' ...
    'ismember("Area1",{''HDB'',''SI'',''VP''})']);  % identified (n = 12)
ChAT = [ChAT 'n067_141017a_1.3'];   % without miss-triggering it passes the cluster criteria
ChAT = [ChAT 'n067_141019a_5.2'];   % light spike assisted clustering

global DATAPATH
impdir = fullfile([DATAPATH,'HDB\taggingsummary_newdata\reliability2\mat_allcells\']);
dr = dir(impdir);
files = {dr(3:end).name};
NumCells = length(ChAT);
for iC = 1:NumCells
    cellid = ChAT{iC};
    cellidt = regexprep(cellid,'\.','_');
    load(fullfile(impdir,[cellidt '_tagging_reliability.mat']))
    if isequal(cellid,'n078_150110a_3.1')   % cell was only there for the first stim session (0 sikes in the other)
        burst_types = burst_types(2);
        Reliability = Reliability(2);
    end
%     figure
%     plot(burst_types/2,Reliability)
    disp([num2str(iC) '  ' cellid])
%     keyboard
    stim_freq{end+1} = burst_types / 2;
%     burst_types / 2
    Rs(end+1,1:length(Reliability)) = Reliability;
end

keyboard