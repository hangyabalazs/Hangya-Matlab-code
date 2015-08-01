%% load NB

load('D:\Dropbox\analysis\NB\average_performance_newdata\PSY_RT_VARS.mat')
GoPerformance1=GoPerformance(1:54);    % n077: only 2 sessions, 10 cells, none cholinergic
NoGoPerformance1=NoGoPerformance(1:54);
GoReactionTime1=GoReactionTime(1:54);
NoGoReactionTime1=NoGoReactionTime(1:54);

%% load HDB

load('D:\Dropbox\analysis\HDB\average_performance_newdata\PSY_RT_VARS.mat')
GoPerformance2=GoPerformance;
NoGoPerformance2=NoGoPerformance;
GoReactionTime2=GoReactionTime;
NoGoReactionTime2=NoGoReactionTime;

%% merge

GoReactionTime = [GoReactionTime1 GoReactionTime2];
NoGoReactionTime = [NoGoReactionTime1 NoGoReactionTime2];
GoPerformance = [GoPerformance1 GoPerformance2];
NoGoPerformance = [NoGoPerformance1 NoGoPerformance2];

%% reshape

GoPerformance = reshape(GoPerformance,3,length(GoPerformance)/3);
NoGoPerformance = reshape(NoGoPerformance,3,length(NoGoPerformance)/3);
GoReactionTime = reshape(GoReactionTime,3,length(GoReactionTime)/3);
NoGoReactionTime = reshape(NoGoReactionTime,3,length(NoGoReactionTime)/3);
SI = [0 0.5 1];

%% plot parameters

gocolor = [0 0.8 0];   % color for go performace
nogocolor = [0.8 0 0];   % color for no-go performance
rtcolor = 'k';   % color for reaction time
ms = 8;   % marker size
lw = 2;   % line width

%% plot

H1 = figure;
set(gcf,'DefaultAxesTickDir','out')
set(gcf,'DefaultAxesBox','off')
xSI = linspace(0,1,length(SI));
plot(SI,mean(GoPerformance,2),'LineStyle','-','Marker','o','MarkerEdgeColor','none',...
    'MarkerFaceColor',gocolor,'Color',gocolor,'MarkerSize',ms,'LineWidth',lw);
hold on
errorbar(SI,mean(GoPerformance,2),nanse(GoPerformance'),'Color',gocolor,'LineWidth',lw)
plot(SI,mean(NoGoPerformance,2),'LineStyle','-','Marker','o','MarkerEdgeColor','none',...
    'MarkerFaceColor',nogocolor,'Color',nogocolor,'MarkerSize',ms,'LineWidth',lw);
errorbar(SI,mean(NoGoPerformance,2),nanse(NoGoPerformance'),'Color',nogocolor,'LineWidth',lw)
ylim([0 1])
box off
title('Performance')
        
H2 = figure;   % reaction time
set(gcf,'DefaultAxesTickDir','out')
set(gcf,'DefaultAxesBox','off')
plot(SI,mean(GoReactionTime,2),'LineStyle','-','Marker','o','MarkerEdgeColor','none',...
    'MarkerFaceColor',rtcolor,'Color',rtcolor,'MarkerSize',ms,'LineWidth',lw);
errorbar(SI,mean(GoReactionTime,2),nanse(GoReactionTime'),'Color',rtcolor,'LineWidth',lw)
box off
title('Reaction time')
xlabel('Sound Pressure Level (dB)')

%% overlay individual curves for each mouse

figure(H1)
hold on
plot(SI,NoGoPerformance,'Color',[1 0.7 0.7])
plot(SI,GoPerformance,'Color',[0.7 1 0.7])
figure(H2)
hold on
plot(SI,GoReactionTime,'Color',[0.7 0.7 0.7])