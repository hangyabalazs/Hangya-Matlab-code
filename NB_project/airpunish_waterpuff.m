function airpunish_waterpuff
%AIRPUNISH_WATERPUFF   False alarm rate.
%   AIRPUNISH_WATERPUFF compares false alarm rate to water+air puff (above
%   0.5 hit rate engagement criterium).
%
%   Se also AIRPUNISH2, AIRPUNISH2B and AIRPUFF_NONOISE.

% Sessions
sessions = {'n081' '150507b'; ...
    'n082' '150507e'; ...
    'n085' '150512a'; ...
    'n086' '150511a'; ...
    'n087' '150510a'};

[ehitrate1, efarate1, ln1] = main(sessions)   % first session (w engaged period >=70; or >50, it's the same - 0 and 5 were skipped)

% Plot
figure;
MN = median([ehitrate1; efarate1],2);
SE = [se_of_median(ehitrate1) se_of_median(efarate1)];
bar(1:2,MN,'BarWidth',0.7,'FaceColor','none',...
    'EdgeColor',[0.8 0 0])
hold on
E = errorbar(1:2,MN,SE,'r+');
errorbar_tick(E,0)

keyboard

% -------------------------------------------------------------------------
function [ehitrate, efarate, ln] = main(sessions)

% False alarm rate
NumSessions = size(sessions,1);
[ehitrate, efarate, ln] = deal(nan(1,NumSessions));
[cumhits cummiss cumfas cumcrs] = deal(0);
for iS = 1:NumSessions
    try
        fnm = fullfile('F:\NB_cellbase\',sessions{iS,1},sessions{iS,2},'TE.mat');
        TE = load(fnm);
    catch
        fnm = fullfile('G:\HDB_cellbase\',sessions{iS,1},sessions{iS,2},'TE.mat');
        TE = load(fnm);
    end
    
    NumTrials = length(TE.TrialStart);
    wn = 20;   % window size (number of trials)
    [hits,fas,gos,nogos,hitrate,farate,restartrate] = deal(nan(NumTrials-1,1));
    for iT = 1:NumTrials-1     % last trial may be incomplete
        
        % Outcome
        if iT >= wn
            hits(iT) = sum(~isnan(TE.Hit(iT-wn+1:iT)));   % number of hits
            fas(iT) = sum(~isnan(TE.FalseAlarm(iT-wn+1:iT)));   % number of false alarms
            gos(iT) = sum(~isnan(TE.Hit(iT-wn+1:iT))|~isnan(TE.Miss(iT-wn+1:iT)));   % number of go tones
            nogos(iT) = sum(~isnan(TE.FalseAlarm(iT-wn+1:iT))|~isnan(TE.CorrectRejection(iT-wn+1:iT)));   % number of no-go tones
            
%             restartrate(iT) = sum(numrestarts(iT-wn+1:iT));   % number of restarts
        end
        
        hitrate(iT) = hits(iT) / gos(iT);   % hit rate
        farate(iT) = fas(iT) / nogos(iT);   % false alarm rate
        
    end
    
    startinx = find(TE.Hit==1,1,'first');   % engagement: from first hit
    endinx = find(hitrate>0.5,1,'last');
    ln(iS) = endinx - startinx;   % period of engagement
    
    ehits = sum(~isnan(TE.Hit(startinx:endinx)));   % number of hits in the engaged period
    efas = sum(~isnan(TE.FalseAlarm(startinx:endinx)));   % number of false alarms in the engaged period
    egos = sum(~isnan(TE.Hit(startinx:endinx))|~isnan(TE.Miss(startinx:endinx)));   % number of go tones in the engaged period
    enogos = sum(~isnan(TE.FalseAlarm(startinx:endinx))|~isnan(TE.CorrectRejection(startinx:endinx)));   % number of no-go tones in the engaged period
    ehitrate(iS) = ehits / egos;   % hit rate
    efarate(iS) = efas / enogos;   % false alarm rate
    
    figure
    plot(hitrate,'g')
    hold on
    plot(farate,'r')
    
    % Chi-square test
    [h p] = b_chi2test2([ehits egos-ehits],[efas enogos-efas],0.05)   % chi square test
    
    % Pooled sample
    endinx = min(endinx,startinx+40);   % keep the sample balanced
    ehits = sum(~isnan(TE.Hit(startinx:endinx)));   % number of hits in the engaged period
    efas = sum(~isnan(TE.FalseAlarm(startinx:endinx)));   % number of false alarms in the engaged period
    egos = sum(~isnan(TE.Hit(startinx:endinx))|~isnan(TE.Miss(startinx:endinx)));   % number of go tones in the engaged period
    enogos = sum(~isnan(TE.FalseAlarm(startinx:endinx))|~isnan(TE.CorrectRejection(startinx:endinx)));   % number of no-go tones in the engaged period
    
    cumhits = cumhits + ehits;
    cummiss = cummiss + egos - ehits;
    cumfas = cumfas + efas;
    cumcrs = cumcrs + enogos - efas;
end
[h p] = b_chi2test2([cumhits cummiss],[cumfas cumcrs],0.05) 
keyboard