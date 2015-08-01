function fig_attentioncells_regression

% Regression between reaction time (dependent variable) and firing rate in
% the foreperiod, but limited to 0.5 s (independent variable), calculated
% by regression_analysis.m
load('C:\Balazs\_analysis\NB\regression\All\new2\regression_results.mat')

% Regression between reaction time (dependent variable) and firing rate in
% the foreperiod, but limited to 0.5 s (independent variable), calculated
% by regression_analysis.m; trials from the lower 15 percentile of the RT
% distribution excluded
% load('C:\Balazs\_analysis\NB\regression\All\impexclude\regression_results.mat')

% Regression between reaction time (dependent variable) and firing rate in
% the foreperiod, but limited to 0.5 s (independent variable), calculated
% by regression_analysis.m; trials from the lower and upper 10 percentile
% of the RT distribution excluded
% load('C:\Balazs\_analysis\NB\regression\All\slowfastexclude\regression_results.mat')


% cellid = 'n023_111214b_7.1';
% cellid = 'n040_121110a_1.2';
% cellid = 'n028_120217a_7.1';  % main fig.2
% cellid = 'n028_120302a_7.1';
cellid = 'n020_111021a_4.1';
regressionplots(cellid)


% -------------------------------------------------------------------------
function regressionplots(cellid,varargin)

% Default arguments
prs = inputParser;
addRequired(prs,'cellid',@iscellid)
g = parse(prs,cellid,varargin{:});

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
[b,stats] = robustfit(x,y);

% Regression plot
figure
plot(x,y,'ko','MarkerfaceColor','k')   % raster plot
xlabel('Firing rate (Hz)')
ylabel('Reaction time (ms)')
setmyplot_Balazs
axis square
icp = b(1);   % intercept
gr = b(2);   % gradient
xx = min(x):0.01:max(x);
yy = xx .* gr + icp;
hold on
plot(xx,yy,'Color',[0.6627 0.6196 0.4039],'LineWidth',2)   % overlay regression line
text('Units','normalized','Position',[0.7 0.7],...
    'String',{['p = ' num2str(p)] ['R = ' num2str(R)]})

% Fit FR on RT - options
s = fitoptions('Method','NonlinearLeastSquares',...   % non-linear least squares
    'Lower',[-Inf -Inf],...  % lower bound
    'Upper',[Inf Inf],...   % upper bound
    'Startpoint',[1 0],...   % initial value
    'Robust','on');  % robust fit

% Fit FR on RT - model
f = fittype('a*x+b','options',s);

% Fit
[fun gof] = fit(x,y,f,s);
figure
plot(y,'k')
hold on
plot(x*fun.a+fun.b,'r')

[fun gof] = fit(y,x,f,s);
figure
plot(x,'k')
hold on
plot(y*fun.a+fun.b,'r')

L5 = plot(fun,'m');
set(L5,'LineWidth',2)


