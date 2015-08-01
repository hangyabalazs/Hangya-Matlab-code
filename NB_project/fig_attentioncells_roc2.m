function fig_attentioncells_roc2
%FIG_ATTENTIONCELLS_ROC2   Attention analysis figures.
%   FIG_ATTENTIONCELLS_ROC2 makes figures for reaction time prediction of
%   cholinergic neurons (ROC window: 50 ms, causal).
%
%   See also REGRESSION_ANALYSIS, STIMRASTER2, NBITIFR and
%   COND_ACCURACY_FR3.

% Directories
global DATAPATH
% inpdir = fullfile(DATAPATH,'NB','attentioncells3_ChAT50ms',filesep);   % results directory
inpdir = fullfile(DATAPATH,'NB','attentioncells3_ChAT_newdata',filesep);   % results directory

% Cell types
ChAT = selectcell(['"ChAT+"==1&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("session_type",{''gonogo'',''feedbackdelay''})&' ...
    'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})']);  % identified ChAT+ cells
ChAT = [ChAT 'n045_121217x_4.6'];   % clustered based on light-evoked spikes

pChAT = selectcell(['"pChAT+"==1&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})']);  % putative ChAT+ cells

cellids = ChAT;
NumCells = length(cellids);   % number of attention cells
ROCs = nan(NumCells,298);
numhits = nan(1,NumCells);
for iC = 1:NumCells   % loop through attention cells
    cellid = cellids{iC};
    cellidt = regexprep(cellid,'\.','_');
    fnm = fullfile(inpdir,[cellidt '_StimulusOn_ROC.mat']);
    
    load(fnm)
    ROCs(iC,:) = smooth(nan2zero(ROC),'linear',5);   % NaN when all spike counts are 0 for both distributions
    
    TE = loadcb(cellid,'TrialEvents');
    numhits(iC) = nansum(TE.Hit);
end
% keyboard

ROCtime = ROCtime + 0.025;   % make the window causal
ROCs2 = ROCs(:,ROCtime>=-2&ROCtime<=0.6);
ROCtime2 = ROCtime(ROCtime>=-2&ROCtime<=0.6);
mn = nanmean(ROCs2,2);
[srt inx] = sort(mn,'descend');   % sort according to mean ROC
figure
imagesc(ROCtime2,1:size(ROCs2,1),ROCs2(inx,:))
set(gca,'CLim',[-0.7 0.7])

figure
% plot(ROCtime2,mean(ROCs2))
errorshade(ROCtime2,nanmean(ROCs2),nanse(ROCs2),'LineColor','k',...
    'ShadeColor','k')
keyboard

% binx = ROCs2(:,end) > 0.1;
% ROCs3 = ROCs2(~binx,:);
% mn = mean(ROCs3,2);
% [srt inx] = sort(mn,'descend');   % sort according to mean ROC
% figure
% imagesc(ROCtime2,1:size(ROCs3,1),ROCs3(inx,:))
% set(gca,'CLim',[-0.7 0.7])
% figure
% plot(ROCtime2,mean(ROCs3))