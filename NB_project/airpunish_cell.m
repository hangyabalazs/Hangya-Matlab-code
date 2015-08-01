%%

% TE=load('Z:\NB_cellbase\n023\111208a\TE.mat');
% TE=load('Z:\NB_cellbase\n029\120131a\TE.mat');
% TE=load('Z:\NB_cellbase\n045\121211a\TE.mat');   % 2nd - min 20 trials of engagement limited by 0.7-crossings
% TE=load('Z:\NB_cellbase\n046\121207a\TE.mat');
% TE=load('Z:\NB_cellbase\n071\141124a\TE.mat');
% TE=load('Z:\NB_cellbase\n072\141124a\TE.mat');

% TE=load('Z:\HDB_cellbase\n067\141001a\TE.mat');  % 3rd
% TE=load('Z:\HDB_cellbase\n070\141108a\TE.mat');
% TE=load('Z:\HDB_cellbase\n078\141207a\TE.mat');

%%

numrestarts = cellfun(@(s)length(s)-1,TE.ITIBegins(1:end-1));
restartinx = cellfun(@(s)length(s)>1&length(s)<7,TE.ITIBegins(1:end-1));   % 1-5 restarts

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
        
        restartrate(iT) = sum(numrestarts(iT-wn+1:iT));   % number of restarts
    end
    
    hitrate(iT) = hits(iT) / gos(iT);   % hit rate
    farate(iT) = fas(iT) / nogos(iT);   % false alarm rate
    
end

%%

figure
plot(hitrate,'go')
hold on
plot(farate,'ro')

%%

figure
plot(TE.Hit,'go','MarkerSize',6,'MarkerFaceColor','g')
hold on
plot(TE.Miss,'ro','MarkerSize',6,'MarkerFaceColor','r')
plot(TE.FalseAlarm-2,'ro','MarkerSize',6,'MarkerFaceColor','r')
plot(TE.CorrectRejection-2,'go','MarkerSize',6,'MarkerFaceColor','g')
ylim([-2 2])

figure;plot(hitrate)
hold on;plot(farate,'r')

