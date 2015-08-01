function incongruent_hits
%INCONGRUENT_HITS   Test for differences in response properties after incongruent hits.
%   INCONGRUENT_HITS compares the number of evoked spikes after congruent
%   and incongruent hits (hits preceded by non-hit trials or hits
%   preceded by false alarms).  Error bars show standard deviations
%   (moderately meaningful).
%
%   See also RT_VS_SPNO.

% Directories
global DATAPATH
resdir = [DATAPATH 'NB\incongruent_hits\'];
issave = true;

% Cholinergic cells
ChAT = selectcell(['"ChAT+"==1&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})']);  % identified
pChAT = selectcell(['"pChAT+"==1&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})']);  % putative
allChAT = [ChAT pChAT];   % identified and putative
NumChAT = length(allChAT);   % number of cholinergic cells

% Number of spikes fired in response to reinforcement
FA_spno = getvalue('FA_SpikeNumberDistribution',allChAT);   % spike numbers for false alarms
Hit_spno = getvalue('Hit_SpikeNumberDistribution',allChAT);   % spike numbers for false alarms

% PSTH statistics
FA_psth_stats = getvalue('FA_psth_stats',allChAT);   % PSTH statistics, false alarms
Hit_psth_stats = getvalue('Hit_psth_stats',allChAT);   % PSTH statistics, hits
FA_psth_stats = nancell2struct(FA_psth_stats);
Hit_psth_stats = nancell2struct(Hit_psth_stats);
FA_act = [FA_psth_stats.Wpa] < 0.01;  % cells activated after false alarms
Hit_act = [Hit_psth_stats.Wpa] < 0.01;   % cells activated after hits
FA_act_inx = find(FA_act);
Hit_act_inx = find(Hit_act);
ChAT_FA_act = allChAT(FA_act);   % cholinergic cells activated by false alarms
ChAT_Hit_act = allChAT(Hit_act);   % cholinergic cells activated by hits
NumFA_act = length(ChAT_FA_act);   % number of cells activated by false alarms
NumHit_act = length(ChAT_Hit_act);   % number of cells activated by hits

% Progress indicator
wb = nicewaitbar(0,'Running ''RT vs fsplatency'' for false alarms...','Name','Please wait...');  % progress indicator
global WB
WB(end+1) = wb;

% Association between previous outcome and cholinergic response
for iC = 1:NumHit_act  % loop through hit-activated cholinergic cells
    inx = Hit_act_inx(iC);
    cellid = allChAT{inx};  % cell ID
    
    % Previous outcome 
    TE = loadcb(cellid,'TrialEvents');  % trial events
    Hits = nan2zero(TE.Hit);   % Hit trials - 1 if the current trial was a hit
    hitinx = find(Hits);
    prevHits = [0 Hits(1:end-1)];   % 1 if prev. trial was a hit
    FAs = nan2zero(TE.FalseAlarm);   % false alarm trials - 1 if the current trial was a false alarm
    prevFAs = [0 FAs(1:end-1)];   % 1 if prev. trial was a false alarm
    incongHits = ~prevHits & Hits;   % incongruent hits: the prev. trial was not a hit
    incongHits2 = prevFAs & Hits;   % the prev. trial was a false alarm
    congHits = prevHits & Hits;   % congruent hits: the prev. trial was also a hit
    gospno = Hit_spno{inx};   % number of spikes fired in response to reinforcement
    
    % Plot association between previous outcome and number of eveoked spikes and return p-value
    figure
    S1 = subplot(1,2,1);
    p = plotassoc(S1,gospno,incongHits(hitinx),congHits(hitinx),'incongruent hits','congruent hits');
    xl = xlim;
    yl = ylim;
    text(xl(1)+diff(xl)*0.7,yl(1)+diff(yl)*0.9,['p-value: ' num2str(p)])
    ylabel('mean spike number')
    S2 = subplot(1,2,2);
    p = plotassoc(S2,gospno,incongHits2(hitinx),congHits(hitinx),'incongruent hits','congruent hits');
    xl = xlim;
    yl = ylim;
    text(xl(1)+diff(xl)*0.7,yl(1)+diff(yl)*0.9,['p-value: ' num2str(p)])
    ylabel('mean spike number')
    if issave
        tt = regexprep(cellid,'\.','_');
        fnm = [resdir tt '_INCONGHITS.fig'];
        saveas(gcf,fnm)   % save scatter plot
    end
    
    close all
    waitbar(iC/NumHit_act)
end
close(wb)

% -------------------------------------------------------------------------
function p = plotassoc(A,X,G1,G2,L1,L2)

% Plot
axes(A);
% boxplot([X(G1) X(G2)],[zeros(1,sum(G1)) ones(1,sum(G2))])
bar([mean(X(G1)) mean(X(G2))],'FaceColor','none','EdgeColor',[0.2 0.8 0.2],'LineWidth',3);
set(A,'XTickLabel',{L1 L2})
hold on
errorbar([1 2],[mean(X(G1)) mean(X(G2))],[std(X(G1)) std(X(G2))],'k+','LineWidth',3)

% Mann-Whitney U-test
p = ranksum(X(G1),X(G2));