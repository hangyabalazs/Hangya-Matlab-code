function depth_vs_faresponse
%DEPTH_VS_FARESPONSE   Correlation between DV position and reward-related firing rate increase.
%   DEPTH_VS_FARESPONSE plots relative firing rate change, i.e. the 
%   logarithm of maximal firing rate divided by baseline firing rate as a 
%   function of dorso-ventral position for all cholinergic and putative
%   cholinergic neurons (see NBCLUSTERING). Correlation coefficient and
%   significance of regression are calculated.
%
%   See also DEPTH_VS_HITRESPONSE, TRAINING_VS_HITRESPONSE, 
%   NBRESPONSESORTER and NBCLUSTERING.

% Cholinergic neurons
ChAT = selectcell(['"ChAT+"==1&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})']);  % identified
ChAT = [ChAT 'n045_121217x_4.6'];   % clustered based on light-evoked spikes
pChAT = selectcell(['"pChAT+"==1&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})']);  % putative
allChAT = [ChAT pChAT];   % identified and putative
NumChAT = length(allChAT);   % number of cholinergic cells

% Color-marker code
mouse = getvalue('RatID',allChAT);  % animal ID
nms = length(unique(mouse));   % number of mice (n=7)
mouse2 = nan(size(mouse));
for k = 1:nms
    mouse2(mouse==min(mouse)) = k;
    mouse(mouse==min(mouse)) = Inf;
end
clr = colormap(jet(nms));   % unique color for each mouse
mrk = 'VSp^hdoX*+<>';   % unique marker for each mouse
colorcode = false;   % control color code
markercode = true;   % control marker code

% Depth
depth = getvalue('DVpos',allChAT);

% Reward related firing rate change
FAstat = getvalue('FA_psth_stats',allChAT);   % hit trial PSTH statistics
FAPSTH = getvalue('FA_psth',allChAT);   % hit trial PSTHs
[FAeff FAeff2 baseline] = deal(nan(1,NumChAT));
figure
hold on
for k = 1:NumChAT
%     disp(allChAT{k})
%     pause(2)
    FAeff(k) = FAstat{k}.maxvalue;   % maximal FR
    baseline(k) = FAstat{k}.baseline;  % baseline FR
    FAeff2(k) = log(FAstat{k}.maxvalue/FAstat{k}.baseline);  % relative FR change on log scale
    if colorcode
        cclr = clr(mouse2(k),:);   % unique color for each mouse
    else
        lc = length(ChAT);   % number of identified neurons
        cclr = [0 0.8*(k<=lc) 0.8*(k>lc)];   % color according to identified/putative
    end
    if markercode
        cmrk = mrk(mouse2(k));   % unique marker for each mouse
    else
        cmrk = 'o';
    end
    plot(depth(k),FAeff2(k),cmrk,'MarkerSize',12,'MarkerFaceColor',cclr,'MarkerEdgeColor','k')  % scatter plot
end

% Regression
y = FAeff2;
x = depth';
[b,bint,r,rint,stats] = regress(y',[ones(length(x),1),x']);
R = sqrt(stats(1));         % correlation coefficient (R-value of the regression)
F = stats(2);           % F-test for H0: all coeff.-s are zero
p = stats(3);           % F-test significance
xlabel('Depth (um)')
ylabel('Response to punishment')
setmyplot_balazs
axis square
icp = b(1);   % intercept
gr = b(2);   % gradient
xx = min(x):0.01:max(x);
yy = xx .* gr + icp;
hold on
plot(xx,yy,'Color','k','LineWidth',2)   % overlay regression line
text('Units','normalized','Position',[0.7 0.7],...
    'String',{['p = ' num2str(p)] ['R = ' num2str(R)]})

keyboard