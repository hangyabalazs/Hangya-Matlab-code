function nbpsthsummary_NT
%NBPSTHSUMMARY   PSTH summary figures.
%   NBPSTHSUMMARY plots population PSTH for activated cholinergic
%   (including putative), activated non-tagged and inhibited non-tagged
%   cells. Z-scored PSTHs are used throughout the analysis. Avergae PSTHs
%   are also plotted after baseline-subtraction.
%
%   Next, average PSTH aligned to reward or punishment is plotted for all
%   (not only significantly modulated) ChAT+ cells. A heat map of the
%   population PSTH for the significantly activated (p<0.05, two-sided
%   Mann-Whitney test) ChAT+ cells is also plotted separately for hit and
%   false alarm trials.
%
%   See also PLOT_POPSTH_CELL_NEW.

%   Balazs Hangya, Cold Spring Harbor Laboratory
%   1 Bungtown Road, Cold Spring Harbor
%   balazs.cshl@gmail.com
%   29-Sept-2013

%   Edit log: BH 9/29/13

% Directories
global DATAPATH
resdir = ([DATAPATH 'NB\psthsummary_NT\']);   % result directory

% Cells
selstr = ['"ChAT+"==1&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("session_type",{''gonogo'',''feedbackdelay''})&' ...
    'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})'];
ChAT = selectcell(selstr);   % cell IDs for ChAT cells
ChAT = [ChAT 'n045_121217x_4.6'];   % clustered based on light-evoked spikes
selstr = ['"pChAT+"==1&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("session_type",{''gonogo'',''feedbackdelay''})&' ...
    'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})'];
pChAT = selectcell(selstr);  % putative
selstr = ['"ChAT+"==0&"pChAT+"==0&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("session_type",{''gonogo'',''feedbackdelay''})&' ...
    'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})'];
NT = selectcell(selstr);   % cell IDs for non-tagged cells (note: some 
% pChATs are included with more than one session here; they can only be
% checked for repetitions AFTER we identify them as pChATs)

% Load data
tags = [ChAT pChAT NT];
allpsth = getvalue('FA_psth',tags);
allpsth = nancell2mat(allpsth);
allpsth = zscore(allpsth,0,2);
allstats = getvalue('FA_psth_stats',tags);
allstats = nancell2struct2(allstats);

% PSTH statistics
activation_peak = [allstats.activation_peak];   % peak time of activation
activation_start = [allstats.activation_start];   % activation onset time
activation_end = [allstats.activation_end];   % activation offset time
maxvalue = [allstats.maxvalue];   % maximal firing rate
Wpa = [allstats.Wpa];   % Mann-Whitney test for significant activation
inhibition_peak = [allstats.inhibition_peak];   % peak time of inhibition
inhibition_start = [allstats.inhibition_start];   % inhibition onset time
inhibition_end = [allstats.inhibition_end];   % inhibition offset time
minvalue = [allstats.minvalue];   % minimal firing rate
Wpi = [allstats.Wpi];   % Mann-Whitney test for significant inhibition
baseline = [allstats.baseline];   % baseline firing rate

% Groups of activated and inhibited cells
inx_act = Wpa < 0.01;   % significant activation
inx_inh = Wpi < 0.01;   % significant inhibition

activated = find(inx_act&~inx_inh);    % indices for activated cells
inhibited = find(inx_inh&~inx_act);    % indices for inhibited cells
ai = find(inx_act&inx_inh);
inx = activation_peak(ai) < inhibition_peak(ai);
activated_inhibited = sort(ai(inx));
inhibited_activated = ai(~inx);
activated = [activated'; activated_inhibited'];   % categorize based on first significant effect
inhibited = [inhibited'; inhibited_activated'];

% Activated NB cells
[~, ~, ChATinx] = intersect(ChAT,tags);   % indices for ChAT cells
[~, ~, pChATinx] = intersect(pChAT,tags);   % indices for pChAT cells
[~, ~, NTinx] = intersect(NT,tags);   % indices for unidentified cells
ChATactinx = intersect(ChATinx,activated);   % indices for activated ChAT cells
pChATactinx = intersect(pChATinx,activated);   % indices for activated pChAT cells
NTactinx = intersect(NTinx,activated);  % indices for activated unidentified cells
NTinhinx = intersect(NTinx,inhibited);  % indices for inhibited unidentified cells
numChAT = length(ChATactinx);   % number of ChAT cells
numNT = length(NTactinx);   % number of unidentified cells
ChAT_PSTH = allpsth(ChATactinx,:);   % PSTHs of ChAT cells
NT_PSTH = allpsth(NTactinx,:);    % PSTHs of NT cells

% [~, ~, ChATinx] = intersect(ChAT,tags);   % indices for ChAT cells
% [~, ~, pChATinx] = intersect(pChAT,tags);   % indices for pChAT cells
% [~, ~, NTinx] = intersect(NT,tags);   % indices for unidentified cells
% allChATinx = [ChATinx pChATinx];   % indices for ChAT and pChAT cells
% ChATactinx = intersect(ChATinx,activated);   % indices for activated ChAT cells
% pChATactinx = intersect(pChATinx,activated);   % indices for activated pChAT cells
% allChATactinx = [ChATactinx pChATactinx];   % indices for activated ChAT and pChAT cells
% NTactinx = intersect(NTinx,activated);  % indices for activated non-ChAT cells
% NTinhinx = intersect(NTinx,inhibited);  % indices for inhibited non-ChAT cells
% NTinhactinx = intersect(NTinx,inhibited_activated);  % indices for inhibited-activated non-ChAT cells

% Average PSTH
wn = [-0.6 0.6];
dt = 0.001;
time = wn(1):dt:wn(2);   % time vector
figure
hold on;
baseline_act = mean(allpsth(NTactinx,1:200));   % baseline-subtracted average z-scores
baseline_inh = mean(allpsth(NTinhinx,1:200));
baseline_pChAT = mean(allpsth(pChATactinx,1:200));
baseline_allNT = nanmean(allpsth(NTinx,1:200));
mn_act = mean(baseline_act);   % baseline for activated NT cells
mn_inh = mean(baseline_inh);   % baseline for inhibited NT cells
mn_pChAT = mean(baseline_pChAT);   % baseline for activated cholinergic cells
mn_allNT = mean(baseline_allNT);   % baseline for all (not only activated) NT cells
errorshade(time,mean(allpsth(NTactinx,:))-mn_act,std(allpsth(NTactinx,:))/sqrt(size(allpsth(NTactinx,:),1)),...
    'LineColor',[0.7 0.7 0.7],'ShadeColor',[0.7 0.7 0.7])
errorshade(time,mean(allpsth(pChATactinx,:))-mn_pChAT,std(allpsth(pChATactinx,:))/sqrt(size(allpsth(pChATactinx,:),1)),...
    'LineColor',[51 204 51]/255,'ShadeColor',[51 204 51]/255)
errorshade(time,mean(allpsth(NTinhinx,:))-mn_inh,std(allpsth(NTinhinx,:))/sqrt(size(allpsth(NTinhinx,:),1)),...
    'LineColor',[1 0.6 0.78],'ShadeColor',[1 0.6 0.78])
errorshade(time,nanmean(allpsth(NTinx,:))-mn_allNT,nanse(allpsth(NTinx,:)),...
    'LineColor',[0.5 0.5 0],'ShadeColor',[0.5 0.5 0])

line([0 0],[-4 20],'Color','k','LineStyle',':')
line([-200 600],[0 0],'Color','k')
set(gca,'box','off',...
    'FontSize',16,'TickDir','out','XLim',[-0.200 0.600],'Ylim',[-2.6 10])
xlabel('Time')
ylabel('Normalized firing rate')
% saveas(gcf,fullfile(resdir,'psth_average.fig'))

% Average PSTHs for identified (not necessarily sign. activated) ChAT cells
green = [51 204 51] / 255;   % colors for plotting
red = [216 41 0] / 255;
NumChAT = length(pChATinx);   % number of identified ChAT+ neurons
FA_psth = getvalue('FA_psth',pChAT);    % PSTH aligned to air puff 
FA_psth = cell2mat(FA_psth);
Hit_psth = getvalue('Hit_psth',pChAT);   % PSTH aligned to water delivery
Hit_psth = cell2mat(Hit_psth);
for k = 1:size(FA_psth,1)   % z-score PSTHs
    FA_psth(k,:) = zscore(FA_psth(k,:));
    Hit_psth(k,:) = zscore(Hit_psth(k,:));
end

inx = time>=-300 & time <= 300;
figure   % plot
errorshade(time,mean(FA_psth),std(FA_psth)/sqrt(size(FA_psth,1)),...
    'LineColor',red,'ShadeColor',red)
hold on
errorshade(time,mean(Hit_psth),std(Hit_psth)/sqrt(size(Hit_psth,1)),...
    'LineColor',green,'ShadeColor',green)
% saveas(gcf,fullfile(resdir,'pChAT_psth_Hit_FA_average.fig'))

[srt Ia] = sort(max(FA_psth,[],2),'descend');   % sort based on FA response
figure   % plot all FA PSTHs, sorted
imagesc(time(inx),1:size(FA_psth,1),FA_psth(Ia,inx))
colormap(hot)
% saveas(gcf,fullfile(resdir,'pChAT_poppsth_FA_all_sorted.fig'))

figure   % plot all Hit PSTHs; sort based on FA response
imagesc(time(inx),1:size(Hit_psth,1),Hit_psth(Ia,inx))
colormap(hot)
% saveas(gcf,fullfile(resdir,'pChAT_poppsth_Hit_all_sorted.fig'))

[m1 m2] = max(FA_psth,[],2);
[srt Ia] = sort(m2,'ascend');   % sort based on FA response
figure   % plot all FA PSTHs, sorted
imagesc(time(inx),1:size(FA_psth,1),FA_psth(Ia,inx))
colormap(hot)
saveas(gcf,fullfile(resdir,'pChAT_poppsth_FA_all_sorted_new.fig'))

figure   % plot all Hit PSTHs; sort based on FA response
imagesc(time(inx),1:size(Hit_psth,1),Hit_psth(Ia,inx))
colormap(hot)
saveas(gcf,fullfile(resdir,'pChAT_poppsth_Hit_all_sorted_new.fig'))

keyboard