function nbattentioncells
%NBATTENTIONCELLS   Attention analysis.
%   NBATTENTIONCELLS performs various analyses to asess whether cells are
%   related to accuracy, action (animal's response), foreperiod
%   distribution or reaction time. See details about the analyses in
%   REGRESSION_ANALYSIS, STIMRASTER2, NBITIFR and COND_ACCURACY_FR3.
%
%   See also REGRESSION_ANALYSIS, STIMRASTER2, NBITIFR and
%   COND_ACCURACY_FR3.

%   Balazs Hangya, Cold Spring Harbor Laboratory
%   1 Bungtown Road, Cold Spring Harbor
%   balazs.cshl@gmail.com
%   30-Oct-2013

%   Edit log: BH 10/30/13

% Regression between reaction time (dependent variable) and firing rate in
% the foreperiod, but limited to 0.5 s (independent variable), calculated
% by regression_analysis.m
% load('C:\Balazs\_analysis\NB\regression\All\new2\regression_results.mat')

% Regression between reaction time (dependent variable) and firing rate in
% the foreperiod, but limited to 0.5 s (independent variable), calculated
% by regression_analysis.m; trials from the lower 15 percentile of the RT
% distribution excluded
% load('C:\Balazs\_analysis\NB\regression\All\impexclude\regression_results.mat')

% Regression between reaction time (dependent variable) and firing rate in
% the foreperiod, but limited to 0.5 s (independent variable), calculated
% by regression_analysis.m; trials from the lower and upper 10 percentile
% of the RT distribution excluded
load('C:\Balazs\_analysis\NB\regression\All\slowfastexclude\regression_results.mat')

% Directories
global DATAPATH
resdir = fullfile(DATAPATH,'NB','attentioncells',filesep);   % results directory
issave = false;

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

VB = selectcell(['"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("Area1",{''VPM'',''VPL'',''VB''})']);
[~, VBinx] = intersect(I,VB);

CPu = selectcell(['"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("Area1",{''CPu'',''Cpu''})']);
[~, CPuinx] = intersect(I,CPu);

RT = selectcell(['"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("Area1",{''RT''})']);
[~, RTinx] = intersect(I,RT);

ACx = selectcell(['"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("Area1",{''ACx''})']);
[~, ACxinx] = intersect(I,ACx);

% Cell types
ChAT = selectcell(['"ChAT+"==1&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})']);  % identified ChAT+ cells
ChAT = [ChAT 'n045_121217x_4.6'];   % clustered based on light-evoked spikes
[~, ChATinx] = intersect(I,ChAT);

pChAT = selectcell(['"pChAT+"==1&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})']);  % putative ChAT+ cells
[~, pChATinx] = intersect(I,pChAT);

PV = selectcell(['"PV+"==1&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})']);  % identified PV+ cells
[~, PVinx] = intersect(I,PV);

% NB cells with significant correlation (p<0.01)
nninx = find(~isnan(p));
NBbinx = intersect(NBinx,nninx);
NBsinx = intersect(NBinx,find(p<0.01));
CORR = I(NBsinx);

% Percentage of significantly activated cells for different areas
pct.NB = sum(p(NBinx)<0.01) / sum(~isnan(p(NBinx)));
pct.RT = sum(p(RTinx)<0.01) / sum(~isnan(p(RTinx)));
pct.VB = sum(p(VBinx)<0.01) / sum(~isnan(p(VBinx)));
pct.ACx = sum(p(ACxinx)<0.01) / sum(~isnan(p(ACxinx)));
pct.CPu = sum(p(CPuinx)<0.01) / sum(~isnan(p(CPuinx)));

% Percentage of significantly activated cells for different cell types
pct.ChAT = sum(p(ChATinx)<0.01) / sum(~isnan(p(ChATinx)));
pct.pChAT = sum(p(pChATinx)<0.01) / sum(~isnan(p(pChATinx)));
pct.PV = sum(p(PVinx)<0.01) / sum(~isnan(p(PVinx)));

% Save
if issave
    fnm = [resdir 'pct.mat'];
    save(fnm,'p','R','pct')
end

% PSTH, raster plot, ROC analysis
nbstimraster2(I(NBsinx),issave)

% ITI cells
nbitifr(I(NBsinx),issave)   % alignment to trial start

NumCells = length(NBsinx);   % number of attention cells
for iC = 1:NumCells   % loop through attention cells
    cellid = I{NBsinx(iC)};
    
    % Regression plots
    regressionplots(cellid,issave,resdir)
    
    % Accuracy cells, action cells
    [H CondPerf] = cond_accuracy_fr3(cellid,'window',[-0.5 0],'limit2iti',true,'event','StimulusOn','display',true);
    cellidt = regexprep(cellid,'\.','_');
    if issave   % save
        fnm = [resdir cellidt '_CondPerf.mat'];
        save(fnm,'CondPerf')
        
        fnm = [resdir cellidt '_AccuracyCRT.fig'];
        saveas(H.Hca,fnm)
        fnm = [resdir cellidt '_AccuracyCRT.jpg'];
        saveas(H.Hca,fnm)
        
        fnm = [resdir cellidt '_ResponseCRT.fig'];
        saveas(H.Hcr,fnm)
        fnm = [resdir cellidt '_ResponseCRT.jpg'];
        saveas(H.Hcr,fnm)
        
        fnm = [resdir cellidt '_RTAccuracy.fig'];
        saveas(H.Hrta,fnm)
        fnm = [resdir cellidt '_RTAccuracy.jpg'];
        saveas(H.Hrta,fnm)
        
        fnm = [resdir cellidt '_FRAccuracy.fig'];
        saveas(H.Hfra,fnm)
        fnm = [resdir cellidt '_FRAccuracy.jpg'];
        saveas(H.Hfra,fnm)
        
        fnm = [resdir cellidt '_FRResponse.fig'];
        saveas(H.Hfrr,fnm)
        fnm = [resdir cellidt '_FRResponse.jpg'];
        saveas(H.Hfrr,fnm)
    end
    close all
end

% -------------------------------------------------------------------------
function regressionplots(cellid,issave,resdir,varargin)

% Default arguments
prs = inputParser;
addRequired(prs,'cellid',@iscellid)
addRequired(prs,'issave',@(s)islogical(s)|ismember(s,[0 1]))
addRequired(prs,'resdir',@ischar)
g = parse(prs,cellid,issave,resdir,varargin{:});

% Load trial events
try
    TE = loadcb(cellid,'TrialEvents');   % load events
    ST = loadcb(cellid,'EVENTSPIKES');   % load prealigned spikes
catch ME
    disp('There was no behavioral protocol for ths session.')
    error(ME.message)
end

% Checking whether 'DeliverFeedback' event is available
sesstype = getvalue('session_type',cellid);
if isequal(sesstype,{'feedbackdelay'})
    alignevent_fa = 'DeliverFeedback';
    alignevent_hit = 'DeliverFeedback';
else
    alignevent_fa = 'LeftPortIn';
    alignevent_hit = 'LeftWaterValveOn';
end

% Relative spike times
stim_pos = findcellstr(ST.events(:,1),'StimulusOn');   % tone onset
response_pos_fa = findcellstr(ST.events(:,1),alignevent_fa);   % animal's respones, false alarms
response_pos_hit = findcellstr(ST.events(:,1),alignevent_hit);   % animal's respones, hits 
spikes_stimulus = ST.event_stimes{stim_pos};   % spikes relative to tone onset
spikes_response = cell(size(spikes_stimulus));
spikes_response(TE.FalseAlarm==1) = ST.event_stimes{response_pos_fa}(TE.FalseAlarm==1);   % spikes relative to response onset, false alarms
spikes_response(TE.Hit==1) = ST.event_stimes{response_pos_hit}(TE.Hit==1);   % spikes relative to response onset, hits

% Calculate regressors and dependent variables
NumTrials = length(TE.TrialStart);   % number of trials
wn = 20;   % window size (number of trials)
[hits,fas,gos,nogos,hitrate,farate,discrimination,dprime,engagement,...
    correct,incorrect,responded,skipped,...
    gort,nogort,gortint,nogortint,gortloc,nogortloc,loudness,iti,...
    fr_stim_0_05,fr_stim_0_10,fr_stim2resp,fr_resp_0_05,fr_resp_0_10,fr_resp_02_0,...
    fr_stim_025_0, fr_stim_05_0, fr_stim_05_02, fr_stim_10_0,...
    fr_iti fr_iti05 fr_iti10] = deal(nan(NumTrials-1,1));
for iT = 1:NumTrials-1     % last trial may be incomplete
    
    % Outcome
    if iT >= wn
        hits(iT) = sum(~isnan(TE.Hit(iT-wn+1:iT)));   % number of hits
        fas(iT) = sum(~isnan(TE.FalseAlarm(iT-wn+1:iT)));   % number of false alarms
        gos(iT) = sum(~isnan(TE.Hit(iT-wn+1:iT))|~isnan(TE.Miss(iT-wn+1:iT)));   % number of go tones
        nogos(iT) = sum(~isnan(TE.FalseAlarm(iT-wn+1:iT))|~isnan(TE.CorrectRejection(iT-wn+1:iT)));   % number of no-go tones
    end
    
    hitrate(iT) = hits(iT) / gos(iT);   % hit rate
    farate(iT) = fas(iT) / nogos(iT);   % false alarm rate
    discrimination(iT) = hitrate(iT) - farate(iT);   % hit - false alarm
    dprime(iT) = norminv(hitrate(iT)) - norminv(farate(iT));   % d' (SDT measure of discrimnability)
    engagement(iT) = hitrate(iT) + farate(iT);   % hit + false alarm
    
    correct(iT) = isequal(TE.Hit(iT),1) | isequal(TE.CorrectRejection(iT),1);
    incorrect(iT) = isequal(TE.FalseAlarm(iT),1) | isequal(TE.Miss(iT),1);
    responded(iT) = isequal(TE.Hit(iT),1) | isequal(TE.FalseAlarm(iT),1);
    skipped(iT) = isequal(TE.Miss(iT),1) | isequal(TE.CorrectRejection(iT),1);
    
    % Reaction time
    gort(iT) = TE.GoRT(iT);   % reaction time
    nogort(iT) = TE.NoGoRT(iT);   % response time for false alarms
    if iT >= wn
        gortint(iT) = nanmean(TE.GoRT(iT-wn+1:iT));   % reaction time averaged over the window
        nogortint(iT) = nanmean(TE.NoGoRT(iT-wn+1:iT));   % no-go response time averaged over the window
    end
    gortloc(iT) = gort(iT) - gortint(iT);   % current reaction time relative to window average
    nogortloc(iT) = nogort(iT) - nogortint(iT);   % current no-go response time relative to window average 
    
    % Sound intensity
    loudness(iT) = TE.StimulusDuration(iT);   % tone intensity
    
    % Foreperiod (expectancy)
    iti(iT) = TE.ITIDistribution(iT);   % length of foreperiod (s)
        
    % Firing rate 0.5 s after stimulus onset
    lim1 = 0;  % start of time window
    lim2 = 0.5;   % end of time window
    lspikes = spikes_stimulus{iT};   % spikes relative to the event
    lspikes2 = lspikes(lspikes>lim1&lspikes<lim2);   % spikes in the time window
    fr_stim_0_05(iT) = length(lspikes2) / (lim2 - lim1);   % firing rate in the time window
    
    % Firing rate 1 s after stimulus onset
    lim1 = 0;  % start of time window
    lim2 = 1;   % end of time window
    lspikes = spikes_stimulus{iT};   % spikes relative to the event
    lspikes2 = lspikes(lspikes>lim1&lspikes<lim2);   % spikes in the time window
    fr_stim_0_10(iT) = length(lspikes2) / (lim2 - lim1);   % firing rate in the time window
    
    % Firing rate in from stimulus to response
    lim1 = 0;
    lim2 = TE.ReactionTime(iT);
    lspikes = spikes_stimulus{iT};   % spikes relative to the event
    lspikes2 = lspikes(lspikes>lim1&lspikes<lim2);   % spikes in the time window
    fr_stim2resp(iT) = length(lspikes2) / (lim2 - lim1);   % firing rate in the time window
    
    % Firing rate 0.25 s before stimulus onset
    lim1 = -0.25;  % start of time window
    lim2 = 0;   % end of time window
    lspikes = spikes_stimulus{iT};   % spikes relative to the event
    lspikes2 = lspikes(lspikes>lim1&lspikes<lim2);   % spikes in the time window
    fr_stim_025_0(iT) = length(lspikes2) / (lim2 - lim1);   % firing rate in the time window
    
    % Firing rate 0.5 s before stimulus onset
    lim1 = -0.5;  % start of time window
    lim2 = 0;   % end of time window
    lspikes = spikes_stimulus{iT};   % spikes relative to the event
    lspikes2 = lspikes(lspikes>lim1&lspikes<lim2);   % spikes in the time window
    fr_stim_05_0(iT) = length(lspikes2) / (lim2 - lim1);   % firing rate in the time window
    
    % Firing rate 0.5-0.2 s before stimulus onset
    lim1 = -0.5;  % start of time window
    lim2 = -0.2;   % end of time window
    lspikes = spikes_stimulus{iT};   % spikes relative to the event
    lspikes2 = lspikes(lspikes>lim1&lspikes<lim2);   % spikes in the time window
    fr_stim_05_02(iT) = length(lspikes2) / (lim2 - lim1);   % firing rate in the time window
    
    % Firing rate 1 s before stimulus onset
    lim1 = -1;  % start of time window
    lim2 = 0;   % end of time window
    lspikes = spikes_stimulus{iT};   % spikes relative to the event
    lspikes2 = lspikes(lspikes>lim1&lspikes<lim2);   % spikes in the time window
    fr_stim_10_0(iT) = length(lspikes2) / (lim2 - lim1);   % firing rate in the time window
    
    % Firing rate 0.5 s after response onset
    lim1 = 0;  % start of time window
    lim2 = 0.5;   % end of time window
    lspikes = spikes_response{iT};   % spikes relative to the event
    lspikes2 = lspikes(lspikes>lim1&lspikes<lim2);   % spikes in the time window
    fr_resp_0_05(iT) = length(lspikes2) / (lim2 - lim1);   % firing rate in the time window
    
    % Firing rate 1 s after response onset
    lim1 = 0;  % start of time window
    lim2 = 1;   % end of time window
    lspikes = spikes_response{iT};   % spikes relative to the event
    lspikes2 = lspikes(lspikes>lim1&lspikes<lim2);   % spikes in the time window
    fr_resp_0_10(iT) = length(lspikes2) / (lim2 - lim1);   % firing rate in the time window
    
    % Firing rate 0.2 s before response onset
    lim1 = -0.2;  % start of time window
    lim2 = 0;   % end of time window
    lspikes = spikes_response{iT};   % spikes relative to the event
    lspikes2 = lspikes(lspikes>lim1&lspikes<lim2);   % spikes in the time window
    fr_resp_02_0(iT) = length(lspikes2) / (lim2 - lim1);   % firing rate in the time window
    
    % Firing rate in the foreperiod
    itiwin = TE.ITIDistribution(iT);   % ITI
    lim1 = -itiwin;
    lim2 = 0;
    lspikes = spikes_stimulus{iT};   % spikes relative to the event
    lspikes2 = lspikes(lspikes>lim1&lspikes<lim2);   % spikes in the time window
    fr_iti(iT) = length(lspikes2) / (lim2 - lim1);   % firing rate in the time window
    
    % Firing rate in the foreperiod, but max 0.5 s
    itiwin = TE.ITIDistribution(iT);   % ITI
    lim1 = max(-itiwin,-0.5);
    lim2 = 0;
    lspikes = spikes_stimulus{iT};   % spikes relative to the event
    lspikes2 = lspikes(lspikes>lim1&lspikes<lim2);   % spikes in the time window
    fr_iti05(iT) = length(lspikes2) / (lim2 - lim1);   % firing rate in the time window
    
    % Firing rate in the foreperiod, but max 1 s
    itiwin = TE.ITIDistribution(iT);   % ITI
    lim1 = max(-itiwin,-1);
    lim2 = 0;
    lspikes = spikes_stimulus{iT};   % spikes relative to the event
    lspikes2 = lspikes(lspikes>lim1&lspikes<lim2);   % spikes in the time window
    fr_iti10(iT) = length(lspikes2) / (lim2 - lim1);   % firing rate in the time window
end

% Regression
% hitinx = intersect(find(~isnan(TE.Hit)),find(gort>prctile(gort,15)&gort<0.6));
hitinx = find(~isnan(TE.Hit));
y = gort(hitinx) * 1000;   % go reaction time
x = fr_iti05(hitinx);   % FR in the foreperiod
X = [ones(length(hitinx),1) x];
[b,bint,r,rint,stats] = regress(y,X);
p = stats(3);
pR = corrcoef(y,X(:,2));
R = pR(3);

% Regression plot
Hr = figure;
plot(x,y,'ko','MarkerfaceColor','k')   % raster plot
xlabel('Firing rate (Hz)')
ylabel('Reaction time (ms)')
setmyplot_balazs
axis square
icp = b(1);   % intercept
gr = b(2);   % gradient
xx = min(x):0.01:max(x);
yy = xx .* gr + icp;
hold on
plot(xx,yy,'Color',[0.6627 0.6196 0.4039],'LineWidth',2)   % overlay regression line
text('Units','normalized','Position',[0.7 0.7],...
    'String',{['p = ' num2str(p)] ['R = ' num2str(R)]})

% Linear prediction of reaction time
Hrtp = figure;
plot(y,'Color',[0.6627 0.6196 0.4039],'LineWidth',2)
hold on
py = X * b;   % linear prediction
plot(py,'Color',[0.2 0.8 0.2],'LineWidth',2)
setmyplot_balazs
legend({'RT' 'linear prediction from FR'})

% Reaction time and inverted firing rate
Hrtfr = figure;
plot(zscore(y),'Color',[0.6627 0.6196 0.4039],'LineWidth',2)
hold on
plot(-zscore(x),'Color',[0.2 0.8 0.2],'LineWidth',2)
setmyplot_balazs
legend({'norm. RT' 'inverted norm. FR'})

% Save
cellidt = regexprep(cellid,'\.','_');
if issave   % save
    fnm = [resdir cellidt '_Regresseion.fig'];
    saveas(Hr,fnm)
    fnm = [resdir cellidt '_Regresseion.jpg'];
    saveas(Hr,fnm)
    
    fnm = [resdir cellidt '_RTPrediction.fig'];
    saveas(Hrtp,fnm)
    fnm = [resdir cellidt '_RTPrediction.jpg'];
    saveas(Hrtp,fnm)
    
    fnm = [resdir cellidt '_RTFR.fig'];
    saveas(Hrtfr,fnm)
    fnm = [resdir cellidt '_RTFR.jpg'];
    saveas(Hrtfr,fnm)
end