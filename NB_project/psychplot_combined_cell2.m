%% load NB

load('D:\Dropbox\analysis\NB\average_performance_newdata\PSY_RT_VARS.mat')
GoPerformance1=GoPerformance;
NoGoPerformance1=NoGoPerformance;
GoReactionTime1=GoReactionTime;
NoGoReactionTime1=NoGoReactionTime;

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

%% overlay individual curves for each mouse

H1 = figure;
hold on
inx = [1:10 12:18 20:22];    % n038 and n077: only 4 and 10 neurons, indiv. plots not shown
plot(SI,NoGoPerformance(:,inx),'Color',[1 0.7 0.7])
plot(SI,GoPerformance(:,inx),'Color',[0.7 1 0.7])

H2 = figure;
hold on
plot(SI,GoReactionTime(:,inx),'Color',[0.7 0.7 0.7])

%% plot

figure(H1)
set(gcf,'DefaultAxesTickDir','out')
set(gcf,'DefaultAxesBox','off')
xSI = linspace(0,1,length(SI));
plot(SI,mean(GoPerformance,2),'LineStyle','-','Color',gocolor,'LineWidth',lw);
hold on
E = errorbar(SI,mean(GoPerformance,2),nanse(GoPerformance'),'Color',gocolor,'LineWidth',lw);
errorbar_tick(E,0);
plot(SI,mean(NoGoPerformance,2),'LineStyle','-','Color',nogocolor,'LineWidth',lw);
E = errorbar(SI,mean(NoGoPerformance,2),nanse(NoGoPerformance'),'Color',nogocolor,'LineWidth',lw);
errorbar_tick(E,0);
ylim([0 1])
box off
title('Performance')
        
figure(H2)   % reaction time
set(gcf,'DefaultAxesTickDir','out')
set(gcf,'DefaultAxesBox','off')
plot(SI,mean(GoReactionTime,2),'LineStyle','-','Color',rtcolor,'LineWidth',lw);
E = errorbar(SI,mean(GoReactionTime,2),nanse(GoReactionTime'),'Color',rtcolor,'LineWidth',lw);
errorbar_tick(E,0);
box off
title('Reaction time')
xlabel('Sound Pressure Level (dB)')