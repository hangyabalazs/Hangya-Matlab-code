function [fa1 fa2 ln] = airpunish

sessions = {'n023' '111208a'; ...
    'n029' '120131a'; ...
    'n045' '121211a'; ...     % 2nd - min 20 trials of engagement limited by 0.7-crossings
    'n046' '121207a'; ...
    'n071' '141125a'; ...   % 2nd
    'n072' '141126a'; ...   % 3rd
    'n067' '141001a'; ...   % 3rd
    'n070' '141108a'; ...
    'n078' '141208a'};   % 2nd
[fa1, fa2, ln] = deal(nan(1,9));
for iS = 1:9
    try
        fnm = fullfile('Z:\NB_cellbase\',sessions{iS,1},sessions{iS,2},'TE.mat');
        TE = load(fnm);
    catch
        fnm = fullfile('Z:\HDB_cellbase\',sessions{iS,1},sessions{iS,2},'TE.mat');
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
    
    fa1(iS) = farate(startinx);   % FA rate, beginning of session
    fa2(iS) = farate(endinx);     % FA rate, end of session
    1;
    
%     figure
%     plot(hitrate,'go')
%     hold on
%     plot(farate,'ro')
%     
%     figure
%     plot(TE.Hit,'go','MarkerSize',6,'MarkerFaceColor','g')
%     hold on
%     plot(TE.Miss,'ro','MarkerSize',6,'MarkerFaceColor','r')
%     plot(TE.FalseAlarm-2,'ro','MarkerSize',6,'MarkerFaceColor','r')
%     plot(TE.CorrectRejection-2,'go','MarkerSize',6,'MarkerFaceColor','g')
%     ylim([-2 2])
%     
%     figure
%     plot(hitrate,'g')
%     hold on
%     plot(farate,'r')
    
end
