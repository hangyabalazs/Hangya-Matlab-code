function RT_vs_fsplatency
%RT_VS_FSPLATENCY   Association between reaction time or stimulus intensity and first spike latency
%   RT_VS_FSPLATENCY calculates regresseion between reaction time or stimulus
%   intensity and latency of the first spike evoked by water or air puff in
%   cholinergic neurons (identified or putative).
%
%   See also RT_VS_SPNO.

% Directories
global DATAPATH
resdir = [DATAPATH 'NB\RT_vs_fsplatency_newdata\'];
issave = true;

% Cholinergic cells
ChAT = selectcell(['"ChAT+"==1&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("session_type",{''gonogo'',''feedbackdelay''})&' ...
    'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})']);  % identified
ChAT = [ChAT 'n045_121217x_4.6'];   % clustered based on light-evoked spikes
pChAT = selectcell(['"pChAT+"==1&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("session_type",{''gonogo'',''feedbackdelay''})&' ...
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

% Association between RT/stimulus intensity and cholinergic response, False Alarms
ngallp = nan(NumFA_act,1);
for iC = 1:NumFA_act  % loop through FA-activated cholinergic cells
    inx = FA_act_inx(iC);
    cellid = allChAT{inx};  % cell ID
    
    % Checking whether 'DeliverFeedback' event is available
    sesstype = getvalue('session_type',cellid);
    if isequal(sesstype,{'feedbackdelay'})
        alignevent = 'DeliverFeedback';
    else
       alignevent = 'LeftPortIn';
    end
    
    % First spike latency
    trialfilter = 'FalseAlarm==1';
    [reliability latency jitter B M lim1 lim2 spikenumberdistribution H] = ...
        reliability_latency_jitter(cellid,...
        'event_type','trial','event',alignevent,'window',[-0.02 0.1],...
        'event_filter','custom','filterinput',trialfilter,'isadaptive',2,...
        'baselinewin',[-0.02 0],'testwin',[0 0.1],'relative_threshold',0.05,...
        'jitterdefinition','burst','display',false);
    valid_trials = filterTrials(cellid,'event_type','trial','event',alignevent,...
        'event_filter','custom','filterinput',trialfilter);  % filter trials
    tsegs_evoked = rel2abstimes(cellid,[lim1 lim2],'trial',alignevent,'valid_trials',valid_trials);   % convert period to epochs relative to event
    NumTrials = length(valid_trials);   % number of valid trials
    fsplatency = nan(1,NumTrials);
    for iT = 1:NumTrials
        selts_evoked = extractSegSpikes(cellid,tsegs_evoked(:,iT));   % find evoked spikes
        if ~isempty(selts_evoked)
            fsplatency(iT) = abs2reltimes(cellid,selts_evoked(1),'trial',alignevent);   % first spike latency
        end
    end
    
    % Reaction time, sound intensity, number of evoked spikes
    TE = loadcb(cellid,'TrialEvents');  % trial events
    ngRT = TE.NoGoRT;   % no-go reaction time
    ngRT = ngRT(~isnan(ngRT));
    ngSI = TE.SoundIntensity;   % stim. intensity
    ngSI = ngSI(~isnan(TE.FalseAlarm));
    ngspno = FA_spno{inx};   % number of spikes fired in response to reinforcement
    
    % Plot association between stim. int. and 1st spike latency and return p-value
    p = plotassoc(ngSI,fsplatency);
    xl = xlim;
    yl = ylim;
    text(xl(1)+diff(xl)*0.7,yl(1)+diff(yl)*0.9,['p-value: ' num2str(p)])
    ylabel('1st spike latency')
    xlabel('Sound intensity')
    if issave
        tt = regexprep(cellid,'\.','_');
        fnm = [resdir tt '_FSPL_VS_NGSI'];
        saveas(gcf,fnm)   % save
    end
    
    % Plot association between RT and 1st spike latency and return p-value
    p = plotassoc2(fsplatency,ngRT);
    xl = xlim;
    yl = ylim;
    text(xl(1)+diff(xl)*0.7,yl(1)+diff(yl)*0.9,['p-value: ' num2str(p)])
    xlabel('1st spike latency')
    ylabel('no-go RT')
    if issave
        tt = regexprep(cellid,'\.','_');
        fnm = [resdir tt '_FSPL_VS_NGRT'];
        saveas(gcf,fnm)   % save scatter plot
    end
    
    % Plot association between number of evoked spikes and 1st spike latency and return p-value
    p = plotassoc(ngspno,fsplatency);
    xl = xlim;
    yl = ylim;
    text(xl(1)+diff(xl)*0.7,yl(1)+diff(yl)*0.9,['p-value: ' num2str(p)])
    ylabel('1st spike latency')
    xlabel('Number of spikes in evoked burst')
    if issave
        tt = regexprep(cellid,'\.','_');
        fnm = [resdir tt '_FSPL_VS_NGSPNO'];
        saveas(gcf,fnm)   % save
    end
    ngallp(inx) = p;  % p-values for regression
    close all
    waitbar(iC/NumFA_act)
end
close(wb)

% Progress indicator
wb = nicewaitbar(0,'Running ''RT vs fsplatency'' for hits...','Name','Please wait...');  % progress indicator
WB(end+1) = wb;

% Association between RT/stimulus intensity and cholinergic response, Hits
goallp = nan(NumHit_act,1);
for iC = 1:NumHit_act  % loop through Hit-activated cholinergic cells
    inx = Hit_act_inx(iC);
    cellid = allChAT{inx};  % cell ID
    
    % Checking whether 'DeliverFeedback' event is available
    sesstype = getvalue('session_type',cellid);
    if isequal(sesstype,{'feedbackdelay'})
        alignevent = 'DeliverFeedback';
    else
       alignevent = 'LeftWaterValveOn';
    end
    
    % First spike latency
    trialfilter = 'Hit==1';
    [reliability latency jitter B M lim1 lim2 spikenumberdistribution H] = ...
        reliability_latency_jitter(cellid,...
        'event_type','trial','event',alignevent,'window',[-0.02 0.1],...
        'event_filter','custom','filterinput',trialfilter,'isadaptive',2,...
        'baselinewin',[-0.02 0],'testwin',[0 0.1],'relative_threshold',0.05,...
        'jitterdefinition','burst','display',false);
    valid_trials = filterTrials(cellid,'event_type','trial','event',alignevent,...
        'event_filter','custom','filterinput',trialfilter);  % filter trials
    tsegs_evoked = rel2abstimes(cellid,[lim1 lim2],'trial',alignevent,'valid_trials',valid_trials);   % convert period to epochs relative to event
    NumTrials = length(valid_trials);   % number of valid trials
    fsplatency = nan(1,NumTrials);
    for iT = 1:NumTrials
        selts_evoked = extractSegSpikes(cellid,tsegs_evoked(:,iT));   % find evoked spikes
        if ~isempty(selts_evoked)
            fsplatency(iT) = abs2reltimes(cellid,selts_evoked(1),'trial',alignevent);   % first spike latency
        end
    end
    
    % Reaction time, sound intensity, number of evoked spikes
    TE = loadcb(cellid,'TrialEvents');  % trial events
    goRT = TE.GoRT;   % go reaction time
    goRT = goRT(~isnan(goRT));
    goSI = TE.SoundIntensity;   % stim. intensity
    goSI = goSI(~isnan(TE.Hit));
    gospno = Hit_spno{inx};   % number of spikes fired in response to reinforcement
    
    % Plot association between stim. int. and 1st spike latency and return p-value
    p = plotassoc(goSI,fsplatency);
    xl = xlim;
    yl = ylim;
    text(xl(1)+diff(xl)*0.7,yl(1)+diff(yl)*0.9,['p-value: ' num2str(p)])
    ylabel('1st spike latency')
    xlabel('Sound intensity')
    if issave
        tt = regexprep(cellid,'\.','_');
        fnm = [resdir tt '_FSPL_VS_GOSI'];
        saveas(gcf,fnm)   % save scatter plot
    end
    
    % Plot association between RT and 1st spike latency and return p-value
    p = plotassoc2(fsplatency,goRT);
    xl = xlim;
    yl = ylim;
    text(xl(1)+diff(xl)*0.7,yl(1)+diff(yl)*0.9,['p-value: ' num2str(p)])
    xlabel('1st spike latency')
    ylabel('go RT')
    if issave
        tt = regexprep(cellid,'\.','_');
        fnm = [resdir tt '_FSPL_VS_GORT'];
        saveas(gcf,fnm)   % save
    end
    
    % Plot association between number of evoked spikes and 1st spike latency and return p-value
    p = plotassoc(gospno,fsplatency);
    xl = xlim;
    yl = ylim;
    text(xl(1)+diff(xl)*0.7,yl(1)+diff(yl)*0.9,['p-value: ' num2str(p)])
    ylabel('1st spike latency')
    xlabel('Number of spikes in evoked burst')
    if issave
        tt = regexprep(cellid,'\.','_');
        fnm = [resdir tt '_FSPL_VS_GOSPNO'];
        saveas(gcf,fnm)   % save
    end
    goallp(inx) = p;  % p-values for regression
    close all
    waitbar(iC/NumHit_act)
end
close(wb)

% -------------------------------------------------------------------------
function p = plotassoc(X,Y)

% inx1 = Y==0;   % no response
% inx2 = Y>=1;   % response
%     if sum(inx1) * sum(inx2) > 0   % neither of the index sets are empty
%         boxstat(ngRT(inx1),ngRT(inx2),'no response','response');  % stat., false alarms
%         title([regexprep(cellid,'_',' ') ' FA'])
%     end

% Association between X and Y
unq = sort(unique(X),'ascend');   % different values of X
NumSI = length(unq);   % number of different values
[xk yk] = deal([]);
[ngmn ngse] = deal(nan(1,NumSI));
for k = 1:NumSI  % loop through different values of X
    sic = unq(k);   % current value
    ngss = Y(X==sic);   % number of Y values corresponding to the given X
    nm = sum(X==sic);
    xk = [xk ngss];   % Y
    yk = [yk ones(1,nm)*k];   % grouping variable
    ngmn(k) = nanmean(ngss);   % mean
    ngse(k) = nanstd(ngss) / sqrt(sum(~isnan(ngss)));   % SE
end

% Plot
figure
subplot(1,2,1)
boxplot(xk,yk)
subplot(1,2,2)
bar(unq,ngmn,0.5,'FaceColor','none')
hold on
errorbar(unq,ngmn,ngse,'k+')

% Regression
y = Y;
x = X;
[b,bint,r,rint,stats] = regress(y',[ones(length(x),1),x']);
R = sqrt(stats(1));         % correlation coefficient (R-value of the regression)
F = stats(2);           % F-test for H0: all coeff.-s are zero
p = stats(3);           % F-test significance
disp(p)

% if p < 0.05
%     plot(unq,ngmn,'k')
% else
%     plot(unq,ngmn,'Color',[0.7 0.7 0.7])
% end

% keyboard

% -------------------------------------------------------------------------
function p = plotassoc2(X,Y)

% Plot
figure
plot(X,Y,'ko','Markersize',12)

% Regression
y = Y;
x = X;
[b,bint,r,rint,stats] = regress(y',[ones(length(x),1),x']);
R = sqrt(stats(1));         % correlation coefficient (R-value of the regression)
F = stats(2);           % F-test for H0: all coeff.-s are zero
p = stats(3);           % F-test significance
disp(p)