function fig_attentioncells_roc
%FIG_ATTENTIONCELLS_ROC   Attention analysis figures.
%   FIG_ATTENTIONCELLS_ROC makes figures for reaction time prediction of
%   significantly predictive neurons (ROC window: 50 ms, causal).
%
%   See also REGRESSION_ANALYSIS, STIMRASTER2, NBITIFR and
%   COND_ACCURACY_FR3.

global DATAPATH
% load('C:\Balazs\_analysis\NB\regression\All\slowfastexclude\regression_results.mat')
load([DATAPATH 'NB\regression_newdata\All\slowfastexclude\regression_results.mat'])

% Directories
% inpdir = fullfile(DATAPATH,'NB','attentioncells3_50ms',filesep);   % results directory
inpdir = fullfile(DATAPATH,'NB','attentioncells3_newdata',filesep);   % results directory

% p-values
NumCells = length(p);
[p2 R2] = deal(nan(1,NumCells));
for k = 1:NumCells
    if ~isempty(p(k).iti05)
        p2(k) = p(k).iti05;
        R2(k) = R(k).iti05;
    end
end
p = p2;
R = R2;

% Areas
NB = selectcell(['"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})']);
[~, NBinx] = intersect(I,NB);

% NB cells with significant correlation (p<0.01)
nninx = find(~isnan(p));
NBbinx = intersect(NBinx,nninx);
NBinx2 = intersect(NBinx,find(p<0.01));
NBsinx = intersect(NBinx2,find(R<0));

cellids = I(NBsinx);
cellids = setdiff(cellids,{'n046_121229d_3.2' 'n046_130102b_2.1'});   % alternative sessions with very few hits (5 and 13)
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
mn = mean(ROCs2,2);
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