function lightpoppsth
%LIGHTPOPPSTH   Population PSTH of light-aligned PSTH.
%   LIGHTPOPPSTH calculates adaptive PSTH (see DAPSTH) aligned to light
%   pulse onset using all pulses of the most efficient stimulation
%   frequency. PSTHs are normalized by peak value.
%
%   See also ULTIMATE_PSTH and DAPSTH.

% Reliability data
global DATAPATH
impdir = fullfile(DATAPATH,'NB\taggingsummary\reliability2\mat\');

% Results directory
resdir = fullfile(DATAPATH,'NB\taggingsummary\poppsth\');

% ChAT cells
selstr = ['"ChAT+"==1&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})'];
ChAT = selectcell(selstr);   % cell IDs for ChAT cells
ChAT = [ChAT 'n045_121217x_4.6'];   % clustered based on light-evoked spikes
NumCells = length(ChAT);

% PSTH
wn = [-0.01 0.02];   % window
dt = 0.001;   % resolution (1 ms)
time = wn(1):dt:wn(2);   % time vector
[psths hsts psths_norm hsts_norm] = deal(nan(NumCells,length(time)));
for iC = 1:NumCells
    cellid = ChAT{iC};   % current cell
    
    % Light stim. frequencies
    SE = loadcb(cellid,'StimEvents');
    burst_types = sort(unique(SE.BurstNPulse),'ascend');
    disp(burst_types)
    
    % Reliability data
    cellidt = regexprep(cellid,'\.','_');
    fnm = fullfile(impdir,[cellidt '_tagging_reliability.mat']);
    load(fnm)
    [~, m2] = max(Reliability);
    bnp = burst_types(m2);
    disp(bnp);
    
    % PSTH
    [psth, spsth, spsth_se, ~, spt] = ...
        ultimate_psth(cellid,'stim','PulseOn',wn,...
        'dt',dt,'display',true,'sigma',0.001,'parts','all','isadaptive',2,...
        'event_filter','custom','filterinput',['BurstNPulse==' num2str(bnp)],'maxtrialno',Inf);   % PSTH
    sspt = sum(spt);
    psths(iC,:) = psth;  % adaptive psth
    hsts(iC,:) = sspt;  % summed biraster
    psths_norm(iC,:) = zscore(psth);  % normalized adaptive psth
    hsts_norm(iC,:) = sspt / max(sspt);  % normalized summed biraster
    
    % Plot
    figure
    bar(time,sspt,'FaceColor','k','EdgeColor','k','BarWidth',1)
    xlim(wn)
end
keyboard

% Plot summary
[m1 m2] = max(hsts_norm,[],2);
[srt Ia] = sort(m2,'ascend');   % sort based on FA response
figure   % plot all FA PSTHs, sorted
imagesc(time,1:NumCells,hsts_norm(Ia,:))
% colormap(hot)
saveas(gcf,fullfile(resdir,'lightpoppsth_sorted.fig'))