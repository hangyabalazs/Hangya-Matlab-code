animalID = 'n023';
animalID2 = 'nb023';
sessions = {'111229a' '111230a' '111231a' '120101a' '120102a' '120103a' ...
    '120104a' '120105a' '120106a'};
numSessions = length(sessions);

PerfGo = cell(1,numSessions);
PerfNoGo = cell(1,numSessions);
TimeGo = cell(1,numSessions);
TimeNoGo = cell(1,numSessions);
Dprime = cell(1,numSessions);
for k = 1:length(sessions)
    sessionID = sessions{k};
    fullpth = [getpref('cellbase','datapath') '\' animalID '\' sessionID '\'];
    TE = load([fullpth 'TE.mat']);

    [PerfGo{k} PerfNoGo{k} TimeGo{k} TimeNoGo{k}] = auditory_gonogo_psychplot2(animalID,sessionID);
    Dprime{k} = PerfGo{k} - PerfNoGo{k};
end

%% mean for nb023

mPerfGo = cell2mat(PerfGo');
mPerfNoGo = cell2mat(PerfNoGo');
mTimeGo = cell2mat(TimeGo');
mTimeNoGo = cell2mat(TimeNoGo');
mDprime = cell2mat(Dprime');

inx = [2:4 7];
ln = length(inx);

PerfGo_nb023 = mean(mPerfGo(inx,1:4));
PerfGoSE_nb023 = std(mPerfGo(inx,1:4)) / sqrt(ln);

PerfNoGo_nb023 = mean(mPerfNoGo(inx,1:4));
PerfNoGoSE_nb023 = std(mPerfNoGo(inx,1:4)) / sqrt(ln);

Dprime_nb023 = mean(mDprime(inx,1:4));
DprimeSE_nb023 = std(mDprime(inx,1:4)) / sqrt(ln);

figure
errorbar([20 30 40 50],PerfGo_nb023,PerfGoSE_nb023,'Color',[0 0.8 0],'LineWidth',3)
hold on
errorbar([20 30 40 50],PerfNoGo_nb023,PerfNoGoSE_nb023,'Color',[0.8 0 0],'LineWidth',3)
set(gca,'LineWidth',2,'Box','off','FontSize',16,'XLim',[15 55])
axis square

figure
errorbar([20 30 40 50],Dprime_nb023,DprimeSE_nb023,'Color',[0 0.8 0],'LineWidth',3)
set(gca,'LineWidth',2,'Box','off','FontSize',16,'XLim',[15 55])
axis square

TimeGo_nb023 = mean(mTimeGo(inx,1:4));
TimeGoSE_nb023 = std(mTimeGo(inx,1:4)) / sqrt(ln);

TimeNoGo_nb023 = mean(mTimeNoGo(inx,1:4));
TimeNoGoSE_nb023 = std(mTimeNoGo(inx,1:4)) / sqrt(ln);

figure
errorbar([20 30 40 50],TimeGo_nb023,TimeGoSE_nb023,'Color',[0 0.8 0],'LineWidth',3)
% hold on
% errorbar([20 30 40 50],TimeNoGo_nb023,TimeNoGoSE_nb023,'Color',[0.8 0 0],'LineWidth',3)
set(gca,'LineWidth',2,'Box','off','FontSize',16,'XLim',[15 55])
axis square