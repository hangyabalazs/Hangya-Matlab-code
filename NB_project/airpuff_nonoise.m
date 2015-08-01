function airpuff_nonoise
%AIRPUFF_NONOISE   False alarm rate.
%   AIRPUFF_NONOISE compares false alarm rate in the first session (above
%   0.7 hit rate engagement criterium: >=70 trials between first and last
%   0.7-crossing) to a later session with comparable difficulty (no
%   background noise).
%
%   Se also AIRPUNISH2 and AIRPUNISH2B.

% Sessions
sessions = {'n074' '141126a'; ...
    'n075' '141126b'; ...   % 2nd
    'n076' '141128a'; ...   % 2nd
    'n079' '141206a'; ...
    'n023' '111208a'; ...
    'n046' '121207a'};

sessions2 = {'n074' '150508a'; ...
    'n075' '150509a'; ...
    'n076' '150508a'; ...
    'n079' '150508a'
    'n023' '111220a'; ...   % step back to 5dB noise
    'n046'  '121229b'};

[ehitrate1, efarate1, ln1] = main(sessions)   % first session (w engaged period >=70; or >50, it's the same - 0 and 5 were skipped)
[ehitrate2, efarate2, ln2] = main(sessions2)   % later session without noise

% Plot
figure;
MN = median([ehitrate1; efarate1; ehitrate2; efarate2],2);
SE = [se_of_median(ehitrate1) se_of_median(efarate1) se_of_median(ehitrate2) se_of_median(efarate2)];
bar(1:4,MN,'BarWidth',0.7,'FaceColor','none',...
    'EdgeColor',[0.8 0 0])
hold on
E = errorbar(1:4,MN,SE,'r+');
errorbar_tick(E,0)

figure;
MN = median([efarate1; efarate2],2);
SE = [se_of_median(efarate1) se_of_median(efarate2)];
bar(1:2,MN,'BarWidth',0.7,'FaceColor','none',...
    'EdgeColor',[0.8 0 0])
hold on
E = errorbar(1:2,MN,SE,'r+');
errorbar_tick(E,0)
axis square
setmyplot_balazs

keyboard

% -------------------------------------------------------------------------
function [ehitrate, efarate, ln] = main(sessions)

% False alarm rate
NumSessions = size(sessions,1);
[ehitrate, efarate, ln] = deal(nan(1,NumSessions));
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
    
    startinx = find(hitrate>=0.7,1,'first');   % engagement: from first to last crossing of 0.7 hit rate
    endinx = find(hitrate>=0.7,1,'last');
    ln(iS) = endinx - startinx;   % period of engagement
    
    ehits = sum(~isnan(TE.Hit(startinx:endinx)));   % number of hits in the engaged period
    efas = sum(~isnan(TE.FalseAlarm(startinx:endinx)));   % number of false alarms in the engaged period
    egos = sum(~isnan(TE.Hit(startinx:endinx))|~isnan(TE.Miss(startinx:endinx)));   % number of go tones in the engaged period
    enogos = sum(~isnan(TE.FalseAlarm(startinx:endinx))|~isnan(TE.CorrectRejection(startinx:endinx)));   % number of no-go tones in the engaged period
    ehitrate(iS) = ehits / egos;   % hit rate
    efarate(iS) = efas / enogos;   % false alarm rate
    
%     figure
%     plot(hitrate,'g')
%     hold on
%     plot(farate,'r')
end