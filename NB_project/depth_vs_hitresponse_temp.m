function depth_vs_hitresponse
%DEPTH_VS_HITRESPONSE   Correlation between DV position and reward-related firing rate increase.
%   DEPTH_VS_HITRESPONSE plots relative firing rate change, i.e. the 
%   logarithm of maximal firing rate divided by baseline firing rate as a 
%   function of dorso-ventral position for all cholinergic and putative
%   cholinergic neurons (see NBCLUSTERING). Correlation coefficient and
%   significance of regression are calculated.
%
%   See also TRAINING_VS_HITRESPONSE, NBRESPONSESORTER and NBCLUSTERING.

% Cholinergic neurons
ChAT = selectcell(['"ChAT+"==1&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})']);  % identified
pChAT = selectcell(['"pChAT+"==1&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})']);  % putative
allChAT = [ChAT pChAT 'n045_121217x_4.6'];   % identified and putative
NumChAT = length(allChAT);   % number of cholinergic cells

% Color code
mouse = getvalue('RatID',allChAT);  % animal ID
nms = length(unique(mouse));   % number of mice (n=7)
mouse2 = nan(size(mouse));
for k = 1:nms
    mouse2(mouse==min(mouse)) = k;
    mouse(mouse==min(mouse)) = Inf;
end
clr = colormap(jet(nms));   % unique color for each mouse

% Depth
depth = getvalue('DVpos',allChAT);

% Reward related firing rate change
Hitstat = getvalue('Hit_psth_stats',allChAT);   % hit trial PSTH statistics
HitPSTH = getvalue('Hit_psth',allChAT);   % hit trial PSTHs
[Hiteff Hiteff2 baseline] = deal(nan(1,NumChAT));
figure
hold on
for k = 1:NumChAT
%     disp(allChAT{k})
%     pause(2)
    Hiteff(k) = Hitstat{k}.maxvalue;   % maximal FR
    baseline(k) = Hitstat{k}.baseline;  % baseline FR
    Hiteff2(k) = log(Hitstat{k}.maxvalue/Hitstat{k}.baseline);  % relative FR change on log scale
    plot(depth(k),Hiteff2(k),'o','MarkerSize',12,'MarkerFaceColor',clr(mouse2(k),:),'MarkerEdgeColor','k')  % scatter plot
end

% Regression
y = Hiteff2;
x = depth';
[b,bint,r,rint,stats] = regress(y',[ones(length(x),1),x']);
R = sqrt(stats(1));         % correlation coefficient (R-value of the regression)
F = stats(2);           % F-test for H0: all coeff.-s are zero
p = stats(3);           % F-test significance

keyboard