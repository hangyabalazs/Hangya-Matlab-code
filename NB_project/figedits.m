%% bimodal ITI

set(gca,'XLim',[0 3],'XTick',[0 1.5 3],'XTickLabel',{'0' '' '3'},'YLim',[0.2 0.4],'YTick',[0.2 0.3 0.4],'YTickLabel',{'200' '' '400'})
xlabel('Foreperiod (s)')
ylabel('Reaction time (ms)')
axis square
setmyplot_balazs

%% bimodal ITI - the distributions

ln = findobj(allchild(gca),'Type','line');
set(ln,'Color',[0.2 0.2 0.2],'LineWidth',1)
set(gca,'XLim',[0 3],'XTick',[0 1.5 3],'XTickLabel',{'1' '' '3'},'YTick',[])
setmyplot_balazs
axis square

%% accuracy cells

axis square
ln = findobj(allchild(gca),'Type','line','LineStyle','-','Color',[0.5 0 0]);
set(ln,'Color',[13 129 64]/255,'LineWidth',1)
ln = findobj(allchild(gca),'Type','line','LineStyle','-','Color',[1 1 0]);
set(ln,'Color',[0 1 0],'LineWidth',1)
ln = findobj(allchild(gca),'Type','line','LineStyle',':','Color',[0.5 0 0]);
set(ln,'LineWidth',1)
ln = findobj(allchild(gca),'Type','line','LineStyle',':','Color',[1 1 0]);
set(ln,'Color',[1 0 0],'LineWidth',1)
legend off
title('')
ylabel('Lick prob.')
set(gca,'XLim',[20 50],'XTick',[20 30 40 50],'YLim',[0 1],'YTick',[0 0.5 1],...
    'YTickLabel',{'0' '' '1'})
setmyplot_balazs

%% accuracy cells - delta-discrimination

tx = findobj(allchild(gca),'Type','text');
delete(tx)
ylabel('\DeltaDiscrimination')
xlabel('SPL (dB)')
title('')
axis square
ln = findobj(allchild(gca),'Type','line');
set(ln,'LineWidth',1)
y_lim = ylim;
set(gca,'XLim',[20 50],'XTick',[20 30 40 50],'YTick',[0 y_lim(2)])
line(xlim,[0 0],'LineStyle',':','Color',[0 0 0],'LineWidth',1)
setmyplot_balazs
set(gcf,'Renderer','painters')

%% RT curves - punishment

cm = hot(6);
clr_new = cm(1:4,:);
clr_old = summer(4);
for k = 1:4
    ln = findobj(allchild(gca),'Type','line','Color',clr_old(k,:));
    set(ln,'LineWidth',1,'Color',clr_new(k,:));
    pt = findobj(allchild(gca),'Type','patch','FaceColor',clr_old(k,:));
    set(pt,'FaceColor',clr_new(k,:));
end
axis([-0.1 0.1 -1 20])
set(gca,'XTick',[-0.1 -0.05 0 0.05 0.1],'XTickLabel',{'-100' '' '0' '' '100'},...
    'YTick',[0 10 20],'YTickLabel',{'0' '' '20'})
legend off
ylabel('Firing rate (Hz)')
xlabel('Time from punsihment (ms)')
setmyplot_balazs
set(gcf,'Renderer','painters')

%% RT curves - reward

cm = hot(6);
clr_new = cm(1:4,:);
clr_old = summer(4);
for k = 1:4
    ln = findobj(allchild(gca),'Type','line','Color',clr_old(k,:));
    set(ln,'LineWidth',1,'Color',clr_new(k,:));
    pt = findobj(allchild(gca),'Type','patch','FaceColor',clr_old(k,:));
    set(pt,'FaceColor',clr_new(k,:));
end
axis([-0.1 0.1 -1 20])
set(gca,'XTick',[-0.1 -0.05 0 0.05 0.1],'XTickLabel',{'-100' '' '0' '' '100'},...
    'YTick',[0 10 20],'YTickLabel',{'0' '' '20'})
legend off
ylabel('Firing rate (Hz)')
xlabel('Time from reward (ms)')
setmyplot_balazs
set(gcf,'Renderer','painters')

%% attention cells - PSTH

xd=get(gco,'xdata');
yd=get(gco,'ydata');
inx=xd<-2|xd>0;
% inx=xd>0;
xd(inx)=[];
yd(inx)=[];
set(gco,'xdata',xd,'ydata',yd)

%% attention cells - PSTH2

axis tight
ln = findobj(allchild(gca),'Type','line','Color',[0 0.5 0]);
set(ln,'LineWidth',1,'Color','k');
ln = findobj(allchild(gca),'Type','line','Color',[0 1 0]);
set(ln,'LineWidth',1,'Color','k','LineStyle',':');
set(gca,'YAxisLocation','right')
pt = findobj(allchild(gca),'Type','patch');
set(pt,'FaceColor','k');
yl = ylim;
ylim([0 round(yl(2)/10)*10])
yl = ylim;
set(gca,'YTick',[yl(1) mean(yl) yl(2)],'YTickLabel',...
    {num2str(yl(1)) '' num2str(yl(2))})
xl = xlim;
set(gca,'XTick',[xl(1) mean(xl) xl(2)],'XTickLabel',...
    {num2str(xl(1)) '' num2str(xl(2))})
xlabel('Time from tone onset (s)')
ylabel({'Firing';'rate (Hz)'})
set(gcf,'Renderer','painters')
setmyplot_balazs

%% attenion cells - regression

ln = findobj(allchild(gca),'Type','line','Color',[0.6627 0.6196 0.4039]);
set(ln,'LineWidth',1,'Color','k');
yl = ylim;
ylim([0 round(yl(2)/100)*100])
yl = ylim;
set(gca,'YTick',[yl(1) mean(yl) yl(2)],'YTickLabel',...
    {num2str(yl(1)) '' num2str(yl(2))})
xl = xlim;
set(gca,'XTick',[xl(1) mean(xl) xl(2)],'XTickLabel',...
    {num2str(xl(1)) '' num2str(xl(2))})
ylabel({'Reaction';'time (ms)'})
setmyplot_balazs

%% raster plots

ln = findobj(allchild(gca),'Type','line','Color','k');
set(ln,'LineWidth',1);
set(gca,'YAxisLocation','left')
set(gca,'Color','w','XColor','k','XTick',[])
line([0 0],ylim,'LineWidth',1,'Color','k')
setmyplot_balazs

%% raster plots - zoomed

ln = findobj(allchild(gca),'Type','line','Color','k');
set(ln,'LineWidth',1);
set(gca,'YAxisLocation','left')
set(gca,'Color','w','XColor','k','XTick',[-0.03 -0.015 0 0.015 0.03],...
    'XTickLabel',{'-30' '' '0' '' '30'},'XLim',[-0.03 0.03])
line([0 0],ylim,'LineWidth',1,'Color','k')
setmyplot_balazs

%% raster plots - zoomed #2

ln = findobj(allchild(gca),'Type','line','Color','k');
set(ln,'LineWidth',1);
set(gca,'YAxisLocation','left')
set(gca,'Color','w','XColor','k','XTick',[-0.05 -0.025 0 0.025 0.05],...
    'XTickLabel',{'-50' '' '0' '' '50'})
line([0 0],ylim,'LineWidth',1,'Color','k')
setmyplot_balazs

%% raster plots - zoomed #3

% figure;viewcell2b('n046_130108a_8.1','TriggerName','DeliverFeedback',...
% 'SortEvent','StimulusOff','eventtype','behav','ShowEvents',...
% {{'StimulusOff'}},'Partitions','#StimulusDuration&Hit','window',[-0.05 0.25],...
% 'Num2Plot',25)
% figure;viewcell2b('n046_130108a_8.1','TriggerName','DeliverFeedback',...
% 'SortEvent','StimulusOff','eventtype','behav','ShowEvents',...
% {{'StimulusOff'}},'Partitions','#StimulusDuration&FalseAlarm','window',[-0.05 0.25],...
% 'Num2Plot',7)

% figure;viewcell2b('n046_130104a_6.2','TriggerName','DeliverFeedback',...
% 'SortEvent','StimulusOff','eventtype','behav','ShowEvents',...
% {{'StimulusOff'}},'Partitions','#StimulusDuration&Hit','window',[-0.05 0.25],...
% 'Num2Plot',24)

% figure;viewcell2b('n046_130104a_6.2','TriggerName','DeliverFeedback',...
% 'SortEvent','StimulusOff','eventtype','behav','ShowEvents',...
% {{'StimulusOff'}},'Partitions','#StimulusDuration&FalseAlarm','window',[-0.05 0.25],...
% 'Num2Plot',15)

ln = findobj(allchild(gca),'Type','line','Color','k');
set(ln,'LineWidth',1);
set(gca,'YAxisLocation','left')
set(gca,'Color','w','XColor','k','XTick',[0 0.05 0.1 0.15 0.2],...
    'XTickLabel',{'0' '' '' '' '200'})
xlim([-0.04 0.24])
line([0 0],ylim,'LineWidth',1,'Color','k')
setmyplot_balazs

%% PSTH

ln = findobj(allchild(gca),'Type','line');
set(ln,'LineWidth',1);
yl = ylim;
ylim([0 round(yl(2)/10)*10])
yl = ylim;
set(gca,'YTick',[yl(1) mean(yl) yl(2)],'YTickLabel',...
    {num2str(yl(1)) '' num2str(yl(2))})
set(gca,'XLim',[-0.3 0.3],'XTick',[-0.3 0 0.3],...
    'XTickLabel',{'-300' '' '300'})
line([0 0],ylim,'LineWidth',1,'Color','k')
legend off
xlabel('Time from reinforcement (ms)')
ylabel('Firing rate (Hz)')
setmyplot_balazs

%% average PSTH

xd=get(gco,'xdata');
yd=get(gco,'ydata');
inx=xd<-300|xd>300;
% inx=xd>0;
xd(inx)=[];
yd(inx)=[];
set(gco,'xdata',xd,'ydata',yd)

%% average PSTH #2

yl = ylim;
set(gca,'XTick',[-300 0 300],'XTickLabel',{'-300' '' '300'},...
    'YTick',[yl(1) 0 yl(2)],'YTickLabel',{num2str(yl(1)) '' num2str(yl(2))})
xlabel('Time from reinforcement (ms)')
ylabel({'Normalized';'firing rate'})
line([0 0],ylim,'LineWidth',1,'Color','k')
setmyplot_balazs
set(gcf,'Renderer','painters')

%% pop PSTH

yl = ylim;
set(gca,'XTick',[-300 0 300],'XTickLabel',{'-300' '' '300'},...
    'YTick',yl(2),'YTickLabel',floor(yl(2)))
xlabel('Time from reinforcement (ms)')
% set(gca,'CLim',[-2 20])
set(gca,'CLim',[-1 20])
setmyplot_balazs

%% pop PSTH zoom in

tmp = get(gco,'cdata');
xd = get(gco,'xdata');
xl = [-300 300];
inx = xd < xl(1) | xd > xl(2);
% inx=xd>0;
tmp(:,inx) = [];
xd(inx) = [];
set(gco,'cdata',tmp,'xdata',xd)
xlim(xl)
set(gca,'XTick',[-300 0 300],'XTickLabel',{'-300' '' '300'})

%% pop PSTH zoom in

tmp = get(gco,'cdata');
xd = get(gco,'xdata');
xl = [-100 100];
inx = xd < xl(1) | xd > xl(2);
% inx=xd>0;
tmp(:,inx) = [];
xd(inx) = [];
set(gco,'cdata',tmp,'xdata',xd)
xlim(xl)
set(gca,'XTick',[-100 0 100],'XTickLabel',{'-100' '' '100'})

%% pop PSTH zoom in

tmp = get(gco,'cdata');
xd = get(gco,'xdata');
xl = [-150 150];
inx = xd < xl(1) | xd > xl(2);
% inx=xd>0;
tmp(:,inx) = [];
xd(inx) = [];
set(gco,'cdata',tmp,'xdata',xd)
xlim(xl)
set(gca,'XTick',[-150 0 150],'XTickLabel',{'-150' '' '150'})

%% box plots

ln = findobj(allchild(gca),'Type','line');
set(ln,'LineWidth',1,'Color','k','LineStyle','-');
om = findobj(allchild(gca),'Tag','Outliers');
set(om,'Marker','none','LineStyle','none')
set(gcf,'Position',[624   457   334   498])
box off
setmyplot_balazs

%% CDF

axis square
set(gca,'YTick',[0 0.5 1],'YTickLabel',{'0' '' '1'})
line(xlim,[0.5 0.5],'Color','black','LineStyle',':','LineWidth',0.75)
setmyplot_balazs

%% Hit vs FA response - latency

axis equal
xlim([0 0.08])
ylim([0 0.04])
box off
ln = findobj(allchild(gca),'Type','line','Color','k');
set(ln,'LineWidth',1,'Color','k','LineStyle',':');
ln = findobj(allchild(gca),'Type','line','MarkerFaceColor',[0 0.8 0]);
set(ln,'MarkerSize',8,'MarkerFaceColor','k','MarkerEdgeColor','k');
set(gca,'XTick',[0 0.02 0.04 0.06 0.08],'XTickLabel',{'0' '' '40' '' '80'},...
    'YTick',[0 0.02 0.04],'YTickLabel',{'0' '20' '40'})
setmyplot_balazs

%% Hit vs FA response - jitter

axis equal
xlim([0 0.01])
ylim([0 0.01])
box off
ln = findobj(allchild(gca),'Type','line','Color','k');
set(ln,'LineWidth',1,'Color','k','LineStyle',':');
ln = findobj(allchild(gca),'Type','line','MarkerFaceColor',[0 0.8 0]);
set(ln,'MarkerSize',8,'MarkerFaceColor','k','MarkerEdgeColor','k');
set(gca,'XTick',[0 0.005 0.01],'XTickLabel',{'0' '' '10'},...
    'YTick',[0 0.005 0.01],'YTickLabel',{'0' '' '10'})
setmyplot_balazs

%% Hit vs FA response - reliability

axis equal
xlim([0 1])
ylim([0 1])
box off
ln = findobj(allchild(gca),'Type','line','Color','k');
set(ln,'LineWidth',1,'Color','k','LineStyle',':');
ln = findobj(allchild(gca),'Type','line','MarkerFaceColor',[0 0.8 0]);
set(ln,'MarkerSize',8,'MarkerFaceColor','k','MarkerEdgeColor','k');
set(gca,'XTick',[0 0.5 1],'XTickLabel',{'0' '' '1'},...
    'YTick',[0 0.5 1],'YTickLabel',{'0' '' '1'})
setmyplot_balazs

%% depth vs FAresponse

% ylim([0 6])
ln = findobj(allchild(gca),'Type','line','MarkerFaceColor',[0 0.8 0]);
set(ln,'MarkerFaceColor',[205 85 160]/255);
ln = findobj(allchild(gca),'Type','line','MarkerFaceColor',[0 0 0.8]);
set(ln,'MarkerFaceColor',[95 215 255]/255);
setmyplot_balazs

%% hit no. vs HITresponse

% ylim([-0.6 4])
xlabel('Number of previous hits')
ylabel('Response to reward')
ln = findobj(allchild(gca),'Type','line','MarkerFaceColor',[0 0.8 0]);
set(ln,'MarkerFaceColor',[205 85 160]/255);
ln = findobj(allchild(gca),'Type','line','MarkerFaceColor',[0 0 0.8]);
set(ln,'MarkerFaceColor',[95 215 255]/255);
% set(gca,'YTick',[0 1 2 3 4])
setmyplot_balazs

%% FA no. vs HITresponse

% ylim([-0.6 4])
xlabel('Number of previous false alarms')
ylabel('Response to reward')
ln = findobj(allchild(gca),'Type','line','MarkerFaceColor',[0 0.8 0]);
set(ln,'MarkerFaceColor',[205 85 160]/255);
ln = findobj(allchild(gca),'Type','line','MarkerFaceColor',[0 0 0.8]);
set(ln,'MarkerFaceColor',[95 215 255]/255);
% set(gca,'YTick',[0 1 2 3 4])
setmyplot_balazs

%% trial no. vs HITresponse

% ylim([-0.6 4])
xlabel('Number of previous trials')
ylabel('Response to reward')
ln = findobj(allchild(gca),'Type','line','MarkerFaceColor',[0 0.8 0]);
set(ln,'MarkerFaceColor',[205 85 160]/255);
ln = findobj(allchild(gca),'Type','line','MarkerFaceColor',[0 0 0.8]);
set(ln,'MarkerFaceColor',[95 215 255]/255);
% set(gca,'YTick',[0 1 2 3 4])
setmyplot_balazs

%% session no. vs HITresponse

% ylim([-0.6 4])
xlabel('Number of previous sessions')
ylabel('Response to reward')
ln = findobj(allchild(gca),'Type','line','MarkerFaceColor',[0 0.8 0]);
set(ln,'MarkerFaceColor',[205 85 160]/255);
ln = findobj(allchild(gca),'Type','line','MarkerFaceColor',[0 0 0.8]);
set(ln,'MarkerFaceColor',[95 215 255]/255);
% set(gca,'YTick',[0 1 2 3 4])
setmyplot_balazs

%% day no. vs HITresponse

% ylim([-0.6 4])
xlabel('Days in training')
ylabel('Response to reward')
ln = findobj(allchild(gca),'Type','line','MarkerFaceColor',[0 0.8 0]);
set(ln,'MarkerFaceColor',[205 85 160]/255);
ln = findobj(allchild(gca),'Type','line','MarkerFaceColor',[0 0 0.8]);
set(ln,'MarkerFaceColor',[95 215 255]/255);
% set(gca,'YTick',[0 1 2 3 4])
setmyplot_balazs

%% light pop PSTH

yl = ylim;
set(gca,'XTick',[-0.01 0 0.01 0.02],'XTickLabel',{'-10' '0' '10' '20'},...
    'YTick',yl(1),'YTickLabel',floor(yl(2)))
xlabel('Time from light onset (ms)')
line([0 0],ylim,'LineWidth',1,'Color','w')
setmyplot_balazs

%% subjective hazard rates

ylim([0 1.2])
box off
ln = findobj(allchild(gca),'Type','line');
set(ln,'Color','k','LineWidth',1);
set(gca,'XTick',[0 1.5 3],'YTick',[0 1.2])
setmyplot_balazs

%% sleep examples - LFP

ln = findobj(allchild(gca),'Type','line');
set(ln,'Color','k','LineWidth',1);
ylim([-10000 8000])
axis off

%% sleep examples - LFP

ln = findobj(allchild(gca),'Type','line');
set(ln,'Color','k','LineWidth',1);
axis off

%% sleep examples - trajectory

ln = findobj(allchild(gca),'Type','line');
set(ln,'Color','k','LineWidth',1);

set(gca,'XTick',[],'YTick',[])
setmyplot_balazs
box on

%% sleep examples - spectrum

ln = findobj(allchild(gca),'Type','line');
set(ln,'Color','k','LineWidth',1);
set(gca,'XTick',[])
setmyplot_balazs

%% tuning curves

xlim([-0.04 0.24])
set(gca,'XTick',[0 0.1 0.2])
ln = findobj(allchild(gca),'Type','line');
set(ln,'LineWidth',2);
legend off
title ''
y_lim = ylim;
set(gca,'XTick',[0 0.1 0.2],'XTickLabel',{'0' '' '200'},'YTick',[0 y_lim(2)/2 y_lim(2)])
setmyplot_balazs

%%

xd=get(gco,'xdata');
yd=get(gco,'ydata');
yds=smooth(yd,'linear',11);
c=get(gco,'color');
hold on;plot(xd,yds,'color',c)

%% HMM model - state

ylim([-0.1 2.1])
ln = findobj(allchild(gca),'Type','line');
set(ln,'LineWidth',1);
set(gca,'XTick',[0 1 2],'YTick',[0 1 2],'YTickLabel',{'start' 'go' 'no-go'})
setmyplot_balazs

%% HMM model - outputs

ylim([0.7 2.3])
ln = findobj(allchild(gca),'Type','line');
set(ln,'Marker','o','MarkerSize',8,'MarkerFaceColor','k')
set(gca,'XTick',[0 1 2],'YTick',[1 2],'YTickLabel',{'G' 'NG'})
setmyplot_balazs

%% HMM model - cumulative outputs

legend off
ylim([0 0.05])
ln = findobj(allchild(gca),'Type','line','Color','g');
set(ln,'LineWidth',1,'Color',[0 0.8 0]);
ln = findobj(allchild(gca),'Type','line','Color','r');
set(ln,'LineWidth',1,'Color',[0.8 0 0]);
set(gca,'XTick',[0 1 2],'YTick',[0 0.05])
setmyplot_balazs

%% HMM model - probabilities

xlim([0 2.38])
ln = findobj(allchild(gca),'Type','line');
set(ln,'LineWidth',1);
set(gca,'XTick',[0 1 2],'YTick',[0 1])
setmyplot_balazs

%% HMM model - licks

xlim([0 2.38])
ylim([0 0.05])
ln = findobj(allchild(gca),'Type','line');
set(ln,'LineWidth',1);
set(gca,'XTick',[0 1 2],'YTick',[0 0.05])
setmyplot_balazs

%% tuning ROC

smk = 7;   % smoothing kernel
load('C:\Balazs\_analysis\NB\tuningcurves_roc_feedback_allChAT_max_25mswindow\n046_130108a_8_1_ROC_fa.mat')
ROCs = smooth(ROC,'linear',smk);
SEs = smooth(SE,'linear',smk);
inx = ROCtime>=-0.04&ROCtime<=0.24;
figure
errorshade(ROCtime(inx),ROCs(inx),SEs(inx),'LineColor',[0.8 0 0],'ShadeColor',[0.8 0 0])
hold on

load('C:\Balazs\_analysis\NB\tuningcurves_roc_feedback_allChAT_max_25mswindow\n046_130108a_8_1_ROC_hit.mat')
ROCs = smooth(ROC,'linear',smk);
SEs = smooth(SE,'linear',smk);
inx = ROCtime>=-0.04&ROCtime<=0.24;
errorshade(ROCtime(inx),ROCs(inx),SEs(inx),'LineColor',[0 0.8 0],'ShadeColor',[0 0.8 0])

xlim([-0.04 0.24])
setmyplot_balazs

%% tuning ROC

load('C:\Balazs\_analysis\NB\tuningcurves_roc_feedback_allChAT_max_25mswindow\n046_130104a_6_1_ROC_fa.mat')
ROCs = smooth(ROC,'linear',smk);
SEs = smooth(SE,'linear',smk);
inx = ROCtime>=-0.04&ROCtime<=0.24;
figure
errorshade(ROCtime(inx),ROCs(inx),SEs(inx),'LineColor',[0.8 0 0],'ShadeColor',[0.8 0 0])
hold on

load('C:\Balazs\_analysis\NB\tuningcurves_roc_feedback_allChAT_max_25mswindow\n046_130104a_6_1_ROC_hit.mat')
ROCs = smooth(ROC,'linear',smk);
SEs = smooth(SE,'linear',smk);
inx = ROCtime>=-0.04&ROCtime<=0.24;
errorshade(ROCtime(inx),ROCs(inx),SEs(inx),'LineColor',[0 0.8 0],'ShadeColor',[0 0.8 0])

xlim([-0.04 0.24])
setmyplot_balazs

%% tuning ROC

load('C:\Balazs\_analysis\NB\tuningcurves_roc_feedback_allChAT_max_25mswindow\n046_130104a_6_2_ROC_fa.mat')
ROCs = smooth(ROC,'linear',smk);
SEs = smooth(SE,'linear',smk);
inx = ROCtime>=-0.04&ROCtime<=0.24;
figure
errorshade(ROCtime(inx),ROCs(inx),SEs(inx),'LineColor',[0.8 0 0],'ShadeColor',[0.8 0 0])
hold on

load('C:\Balazs\_analysis\NB\tuningcurves_roc_feedback_allChAT_max_25mswindow\n046_130104a_6_2_ROC_hit.mat')
ROCs = smooth(ROC,'linear',smk);
SEs = smooth(SE,'linear',smk);
inx = ROCtime>=-0.04&ROCtime<=0.24;
errorshade(ROCtime(inx),ROCs(inx),SEs(inx),'LineColor',[0 0.8 0],'ShadeColor',[0 0.8 0])

%% tuning ROC, newdata

load('d:\Dropbox\analysis\NB\tuningcurves_roc_newdata_feedback_allChAT_max\n046_130102x_4_1_ROC_fa.mat')
ROCs = smooth(ROC,'linear',smk);
SEs = smooth(SE,'linear',smk);
inx = ROCtime>=-0.04&ROCtime<=0.24;
figure
errorshade(ROCtime(inx),ROCs(inx),SEs(inx),'LineColor',[0.8 0 0],'ShadeColor',[0.8 0 0])
hold on

load('d:\Dropbox\analysis\NB\tuningcurves_roc_newdata_feedback_allChAT_max\n046_130102x_4_1_ROC_hit.mat')
ROCs = smooth(ROC,'linear',smk);
SEs = smooth(SE,'linear',smk);
inx = ROCtime>=-0.04&ROCtime<=0.24;
errorshade(ROCtime(inx),ROCs(inx),SEs(inx),'LineColor',[0 0.8 0],'ShadeColor',[0 0.8 0])

xlim([-0.04 0.24])
setmyplot_balazs

%% tuning ROC sign. - newdata

uiopen('d:\Dropbox\analysis\NB\tuningcurves_roc_newdata_feedback_allChAT_max\n046_130102x_4_1_ROCP_hit.fig',1)
xd=get(gco,'xdata');
yd=get(gco,'ydata');
figure;plot(xd,yd<0.05)
xlim([-0.04 0.24])

%% tuning ROC sign.

uiopen('C:\Balazs\_analysis\NB\tuningcurves_roc_feedback_allChAT_max_25mswindow\n046_130108a_8_1_ROCP_hit.fig',1)
xd=get(gco,'xdata');
yd=get(gco,'ydata');
figure;plot(xd,yd<0.01)
xlim([-0.04 0.24])

%% surprise summary bargraph

set(gca,'XTick',[],'YTick',[0 0.5 1])
setmyplot_balazs
set(gco,'FaceColor','none')


%% mean tuning ROC - hit

% load('C:\Balazs\_analysis\NB\tuningcurves_roc_feedback_allChAT_max1\vars.mat','outcellids_hit')
load('d:\Dropbox\analysis\NB\tuningnew_feedback_ChAT_max\vars.mat','outcellids_hit')
cellids = outcellids_hit;
NumCells = length(cellids);
ROCs = nan(NumCells,116);
ROCns = nan(NumCells,116);
Ps = nan(NumCells,116);
numhits = nan(1,NumCells);
% inpdir = 'C:\Balazs\_analysis\NB\tuningcurves_roc_feedback_allChAT_max_25mswindow\';
inpdir = 'd:\Dropbox\analysis\NB\tuningcurves_roc_newdata_feedback_allChAT_max\';
for iC = 1:NumCells   % loop through cells
    cellid = cellids{iC};
    cellidt = regexprep(cellid,'\.','_');
    fnm = fullfile(inpdir,[cellidt '_ROC_hit.mat']);
    
    load(fnm)
    ROCns(iC,:) = ROC;
    ROCs(iC,:) = smooth(nan2zero(ROC),'linear',smk);   % NaN when all spike counts are 0 for both distributions
    Ps(iC,:) = P;
    
%     fnm = fullfile(inpdir,[cellidt '_ROC_hit.fig']);
%     uiopen(fnm,1)
%     fnm = fullfile(inpdir,[cellidt '_ROCP_hit.fig']);
%     uiopen(fnm,1)
    
%     TE = loadcb(cellid,'TrialEvents');
%     numhits(iC) = nansum(TE.Hit);
end
ROC = nanmean(ROCs);
SE = nanse(ROCs);

Hroc = figure;   % plot
inx = ROCtime>=-0.04 & ROCtime<=0.24;
ROCs2 = ROCs(:,inx);
Ps2 = Ps(:,inx);
errorshade(ROCtime(inx),ROC(inx),SE(inx),...
    'LineColor','k','ShadeColor','k')

[jnk, pmn] = nanmax(ROCs2.*zero2nan(double((Ps2<0.05))),[],2);
[jnk, mn] = sort(pmn);
[srt inx2] = sort(mn,'descend');   % sort according to mean ROC

figure
imagesc(ROCtime(inx),1:size(ROC,1),ROCs(mn,inx))
% imagesc(ROCtime(inx),1:size(ROC,1),ROCs(:,inx).*zero2nan(double((Ps(:,inx)<0.05))))
% imagesc(ROCtime(inx),1:size(ROC,1),ROCs2(inx2,:).*zero2nan(double((Ps2(inx2,:)<0.05))))
caxis([-0.25 0.25])

figure
imagesc(ROCtime(inx),1:size(ROC,1),ROCs(1:17,inx))
caxis([-0.25 0.25])

figure
I = ROCs > 0 & Ps < 0.05;
inx2 = ROCtime>0 & ROCtime<0.1;
imagesc(ROCtime(inx2),1:size(ROC,1),I(1:17,inx2))

Pv = nan(1,116);
for k = 1:116
    Pv(k) = signrank(ROCns(:,k));
end
figure;plot(ROCtime,Pv)

%% mean tuning ROC - false alarm

% load('C:\Balazs\_analysis\NB\tuningcurves_roc_feedback_allChAT_max_25mswindow\vars.mat')
cellids = outcellids_hit;
NumCells = length(cellids);   % number of attention cells
ROCs = nan(NumCells,116);
ROCns = nan(NumCells,116);
Ps = nan(NumCells,116);
numhits = nan(1,NumCells);
% inpdir = 'C:\Balazs\_analysis\NB\tuningcurves_roc_feedback_allChAT_max_25mswindow\';
inpdir = 'd:\Dropbox\analysis\NB\tuningcurves_roc_newdata_feedback_allChAT_max\';
for iC = 1:NumCells   % loop through attention cells
    cellid = cellids{iC};
    cellidt = regexprep(cellid,'\.','_');
    fnm = fullfile(inpdir,[cellidt '_ROC_fa.mat']);
    
    load(fnm)
    ROCns(iC,:) = ROC;
    ROCs(iC,:) = smooth(nan2zero(ROC),'linear',smk);   % NaN when all spike counts are 0 for both distributions
    Ps(iC,:) = P;
%     TE = loadcb(cellid,'TrialEvents');
%     numhits(iC) = nansum(TE.Hit);
end
ROC = nanmean(ROCs);
SE = nanse(ROCs);

Hroc = figure;   % plot
inx = ROCtime>=-0.04 & ROCtime<=0.24;
ROCs2 = ROCs(:,inx);
Ps2 = Ps(:,inx);
errorshade(ROCtime(inx),ROC(inx),SE(inx),...
    'LineColor','k','ShadeColor','k')

% mn = nanmean(ROCs2,2);
% [srt inx2] = sort(mn,'descend');   % sort according to mean ROC

% [jnk, pmn] = nanmax(ROCs2.*zero2nan(double((Ps2<0.05))),[],2);
% [jnk, mn] = sort(pmn);

figure
% imagesc(ROCtime(inx),1:size(ROC,1),ROCs(mn,inx))
imagesc(ROCtime(inx),1:size(ROC,1),ROCs(mn,inx))
% imagesc(ROCtime(inx),1:size(ROC,1),ROCs2(inx2,:).*zero2nan(double((Ps2(inx2,:)<0.05))))
caxis([-0.25 0.25])

Pv = nan(1,116);
for k = 1:116
    Pv(k) = signrank(ROCns(:,k));
end
figure;plot(ROCtime,Pv)

%% mean tuning ROC sign.

uiopen('C:\Balazs\_analysis\NB\tuningcurve_summary2\meanROC_hit_pvalue.fig',1)
xd=get(gco,'xdata');
yd=get(gco,'ydata');
figure;plot(xd,yd<0.001)
xlim([-0.04 0.24])   % hit: sign 25-60 ms p<0.001; FA: not sign.

%% reliability, latency, jitter for naive responses

animalID = 'nb053';
sessionID = '140430b';
cellid = [animalID '_' sessionID '_2.2'];   % air puff

alignevent = 'DeliverFeedback';
choosecb('pavlovian')
trialfilter = 'PTrial==1';

[reliability latency jitter B M lim1 lim2 spikenumberdistribution H] = ...
    reliability_latency_jitter(cellid,...
    'event_type','trial','event',alignevent,'window',[-0.02 0.1],...
    'event_filter','custom','filterinput',trialfilter,'isadaptive',2,...
    'baselinewin',[-0.02 0],'testwin',[0 0.1],'relative_threshold',0.05,...
    'jitterdefinition','burst','display',true);

%% reliability, latency, jitter for naive responses

animalID = 'nb053';
sessionID = '140502a';
cellid = [animalID '_' sessionID '_2.4'];   % shock

alignevent = 'DeliverFeedback';
choosecb('pavlovian')
trialfilter = 'PTrial==1';

[reliability latency jitter B M lim1 lim2 spikenumberdistribution H] = ...
    reliability_latency_jitter(cellid,...
    'event_type','trial','event',alignevent,'window',[-0.02 0.1],...
    'event_filter','custom','filterinput',trialfilter,'isadaptive',2,...
    'baselinewin',[-0.02 0],'testwin',[0 0.1],'relative_threshold',0.05,...
    'jitterdefinition','burst','display',true);

%% reliability, latency, jitter for naive responses

animalID = 'h006';
sessionID = '140409c';
cellid = [animalID '_' sessionID '_1.1'];   % air puff

alignevent = 'DeliverFeedback';
choosecb('HDB')
trialfilter = 'PTrial==1';

[reliability latency jitter B M lim1 lim2 spikenumberdistribution H] = ...
    reliability_latency_jitter(cellid,...
    'event_type','trial','event',alignevent,'window',[-0.02 0.1],...
    'event_filter','custom','filterinput',trialfilter,'isadaptive',2,...
    'baselinewin',[-0.02 0],'testwin',[0 0.1],'relative_threshold',0.05,...
    'jitterdefinition','burst','display',true);

%% pop PSTH - HDB

xlim([-0.3 0.3])
axis square
yl = ylim;
set(gca,'XTick',[-300 0 300],'XTickLabel',{'-300' '' '300'},...
    'YTick',yl(2),'YTickLabel',floor(yl(2)))
xlabel('Time from reinforcement (ms)')
set(gca,'CLim',[-2 20])
setmyplot_balazs

%% pChAT clustering

axis square
set(gca,'XTick',[-400 0 400 800],'YTick',[-400 0 400])
setmyplot_balazs

%% pChAT clustering PSTH

xd=get(gco,'xdata');
yd=get(gco,'ydata');
inx=xd<-300|xd>300;
% inx=xd>0;
xd(inx)=[];
yd(inx)=[];
set(gco,'xdata',xd,'ydata',yd)

%% pChAT clustering PSTH #2

yl = ylim;
set(gca,'XTick',[-300 0 300],'XTickLabel',{'-300' '' '300'},...
    'YTick',[yl(1) 0 yl(2)],'YTickLabel',{num2str(yl(1)) '' num2str(yl(2))})
xlabel('Time from reinforcement (ms)')
ylabel({'Normalized';'firing rate'})
line([0 0],ylim,'LineWidth',1,'Color','k')
setmyplot_balazs
set(gcf,'Renderer','painters')

%% blink PSTH

xd=get(gco,'xdata');
yd=get(gco,'ydata');
inx=xd<-0.055|xd>0.055;
% inx=xd>0;
xd(inx)=[];
yd(inx)=[];
set(gco,'xdata',xd,'ydata',yd)

%% RT

ylim([0.14 0.34])
set(gca,'XTick',[0 0.5 1],'YTick',[0.14 0.24 0.34],'YTickLabel',{'140' '' '340'})
axis square
setmyplot_balazs

%% psych performance

set(gca,'XTick',[0 0.5 1],'YTick',[0 0.5 1])
axis square
setmyplot_balazs

%% lick PSTH

xd=get(gco,'xdata');
yd=get(gco,'ydata');
inx=xd<-0.4|xd>1;
% inx=xd>0;
xd(inx)=[];
yd(inx)=[];
set(gco,'xdata',xd,'ydata',yd)

%% lick PETH

xlim([-0.4 1])
axis square
set(gca,'XTick',[0 0.5 1])
setmyplot_balazs

%% tone responses

xd=get(gco,'xdata');
yd=get(gco,'ydata');
yds=smooth(yd,'linear',101);
c=get(gco,'color');
hold on;plot(xd,yds,'color',c)

%% tone responses #2

xlim([-0.4 1])
ylim([0 1])
set(gca,'XTick',[0 0.5 1])
ln = findobj(allchild(gca),'Type','line');
set(ln,'LineWidth',2);
legend off
title ''
y_lim = ylim;
set(gca,'YTick',[0 y_lim(2)/2 y_lim(2)])
setmyplot_balazs