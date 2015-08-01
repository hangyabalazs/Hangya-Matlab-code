function [fa1 fa2 ln] = airpunish2b
%AIRPUNISH2B   False alarm rate.
%   [FA1 FA2 LN] = AIRPUNISH2 compares false alarm rate in the beginning
%   (FA1) and end (FA2) of first session. An engagement driterion is used:
%   the session limits are determined by first and last 0.7-crossings of
%   hit rate. Only sessions where these events span at least 50 trials
%   (LN>=50) are considered when finding first sessions. AIRPUNISH2B
%   analyzes all mice that have sufficient behavioral data (n = 19).
%
%   See also AIRPUNISH2.

% Sessions
sessions = {'n013' '110504a'; ...
    'n018' '110914a'; ...   % 3rd
    'n020' '111010a'; ...   % 2nd
    'n023' '111208a'; ...
    'n026' '120128a'; ...
    'n027' '120128a'; ...
    'n028' '120131a'; ...
    'n029' '120131a'; ...
    'n037' '120830a'; ...
    'n038' '120901a'; ...
    'n039' '120904a'; ...   % 2nd
    'n043' '121125a'; ...
    'n045' '121211a'; ...     % 2nd - min 50 trials of engagement limited by 0.7-crossings
    'n046' '121207a'; ...
    'n071' '141125a'; ...   % 2nd
    'n072' '141126a'; ...   % 3rd
%     'n077' '141129b'; ...   % 2nd
    'n067' '141001a'; ...   % 3rd
    'n070' '141108a'; ...
    'n078' '141208a'};   % 2nd

% False alarm rate
[fa1, fa2, ln] = deal(nan(1,9));
for iS = 1:size(sessions,1)
    disp(iS)
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
    figure
    plot(TE.Hit,'go','MarkerSize',6,'MarkerFaceColor','g')
    hold on
    plot(TE.Miss,'ro','MarkerSize',6,'MarkerFaceColor','r')
    plot(TE.FalseAlarm-2,'ro','MarkerSize',6,'MarkerFaceColor','r')
    plot(TE.CorrectRejection-2,'go','MarkerSize',6,'MarkerFaceColor','g')
    ylim([-2 2])
    
    figure
    plot(hitrate,'g')
    hold on
    plot(farate,'r')
    
end

% Load NB learning curves
load('D:\Dropbox\analysis\NB\learning_curve_newdata\MEAN_LC.mat')  % n077 not saved - only 2 sessions, 10 cells, none cholinergic
GoPerformance1 = aGoPerformance;    % n077: only 2 sessions, 10 cells, none cholinergic
NoGoPerformance1 = aNoGoPerformance;

% Load HDB learning curves
load('D:\Dropbox\analysis\HDB\learning_curve_newdata\MEAN_LC.mat')
GoPerformance2 = aGoPerformance;
NoGoPerformance2 = aNoGoPerformance;

% Merge
aGoPerformance = [GoPerformance1; GoPerformance2];
aNoGoPerformance = [NoGoPerformance1; NoGoPerformance2];

% Means
mn = nanmedian(aNoGoPerformance(:,10:15),2);   % mean FA rate from 10 to 15 sessions
fa3 = mn(~isnan(mn))';   % late FA rate
iinx = [1 3 2 1 1 1 1 1 1 1 2 1 2 1 2 3 3 1 2];
mn = nanmedian(aNoGoPerformance(:,iinx:iinx+4),2);   % mean FA rate from 1 to 5 sessions
fa0 = mn(~isnan(mn))';   % early FA rate

% Plot
figure;
grey = [0.7 0.7 0.7];
bar(1:3,median([fa1; fa2; fa3],2),'BarWidth',0.7,'FaceColor','none',...
    'EdgeColor',[0.8 0 0])
hold on
SE = [se_of_median(fa1) se_of_median(fa2) se_of_median(fa3)];
E = errorbar(1:3,median([fa1; fa2; fa3],2),SE,'r+');
errorbar_tick(E,0)
axis square
setmyplot_balazs

signrank(fa1,fa2)
signrank(fa1,fa3)
signrank(fa2,fa3)

keyboard