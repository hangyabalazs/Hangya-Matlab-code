%% load NB

load('D:\Dropbox\analysis\NB\learning_curve_newdata\MEAN_LC.mat')  % n077 not saved - only 2 sessions, 10 cells, none cholinergic
GoPerformance1 = aGoPerformance;    % n077: only 2 sessions, 10 cells, none cholinergic
NoGoPerformance1 = aNoGoPerformance;

%% load HDB

load('D:\Dropbox\analysis\HDB\learning_curve_newdata\MEAN_LC.mat')
GoPerformance2 = aGoPerformance;
NoGoPerformance2 = aNoGoPerformance;


%% merge

aGoPerformance = [GoPerformance1; GoPerformance2];
aNoGoPerformance = [NoGoPerformance1; NoGoPerformance2];

%% Grand average

mGoPerformance = nanmean(aGoPerformance);   % mean
mNoGoPerformance = nanmean(aNoGoPerformance);
seGoPerformance = nanse(aGoPerformance);   % SEM
seNoGoPerformance = nanse(aNoGoPerformance);

%% plot parameters

gocolor = [0 0.8 0];   % color for go performace
nogocolor = [0.8 0 0];   % color for no-go performance
rtcolor = 'k';   % color for reaction time
ms = 8;   % marker size
lw = 2;   % line width

%% plot

% Plot
H1 = figure;   % performance
set(gcf,'DefaultAxesTickDir','out')
set(gcf,'DefaultAxesBox','off')
plot(mGoPerformance,'LineStyle','-','Marker','o','MarkerEdgeColor','none',...
    'MarkerFaceColor',gocolor,'Color',gocolor,'MarkerSize',ms,'LineWidth',lw);
hold on
E = errorbar(mGoPerformance,seGoPerformance,'Color',gocolor,'LineWidth',lw);
errorbar_tick(E,0)   % eliminate horizontal line from errorbar
plot(mNoGoPerformance,'LineStyle','-','Marker','o','MarkerEdgeColor','none',...
    'MarkerFaceColor',nogocolor,'Color',nogocolor,'MarkerSize',ms,'LineWidth',lw);
E = errorbar(mNoGoPerformance,seNoGoPerformance,'Color',nogocolor,'LineWidth',lw);
errorbar_tick(E,0)   % eliminate horizontal line from errorbar
ylim([0 1])
box off
title('Learning curve')
xlabel('%lick')
ylabel('Session')

%% overlay individual curves for each mouse

figure(H1)
hold on
plot(aGoPerformance','Color',[0.7 1 0.7])
plot(aNoGoPerformance','Color',[1 0.7 0.7])