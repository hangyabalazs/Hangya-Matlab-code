%% load

load('C:\Balazs\_analysis\VIP\psth_variables.mat')

%% groups

tagged =  logical(vldty) & (isact==2);
inx_act = logical(vldty) & (isact==1) & (baseline>2);
inx_inh = logical(vldty) & (isinh) & (baseline>2);
activated = setdiff(find(inx_act&~inx_inh),tagged);
inhibited = find(inx_inh&~inx_act);
ai = find(inx_act&inx_inh);
inx = activation_peak(ai) < inhibition_peak(ai);
activated_inhibited = sort([ai(inx); 10]);    % 10 is activated-inhibited
inhibited_activated = setdiff(ai(~inx),10);
activated = [activated; activated_inhibited];
inhibited = [inhibited; inhibited_activated];

%% peak time

figure
semilogy(activation_peak(activated),maxvalue(activated)./baseline(activated),'o','MarkerEdgeColor',[0.32 0.19 0.19],'MarkerFaceColor',[0.32 0.19 0.19],'MarkerSize',8)
hold on
semilogy(inhibition_peak(inhibited),minvalue(inhibited)./baseline(inhibited),'o','MarkerEdgeColor',[0.48 0.06 0.89],'MarkerFaceColor',[0.48 0.06 0.89],'MarkerSize',8)
semilogy(activation_peak(tagged),maxvalue(tagged)./baseline(tagged),'o','MarkerEdgeColor',[0 0.7 0],'MarkerFaceColor',[0 0.7 0],'MarkerSize',8)
line([0 70],[1 1],'Color','k')
set(gca,'box','off','XTick',0:10:70,'YTick',[0.001 0.01 0.1 1 10 100 1000],'YMinorTick','off',...
    'FontSize',16,'TickDir','out','XLim',[0 71],'Ylim',[0.001 1000])
xlabel('Peak time')
ylabel('Relative change in firing rate')


%% pdf #3 - accepted version

edges = [0 1 3 5 7:3:100];
centers = (edges(1:end-1) + edges(2:end)) / 2;
dist_activated = histc(activation_peak(activated),edges);
dist_activated = dist_activated(1:end-1);
dist_activated = dist_activated / sum(dist_activated);

dist_inhibited = histc(inhibition_peak(inhibited),edges);
dist_inhibited = dist_inhibited(1:end-1);
dist_inhibited = dist_inhibited / sum(dist_inhibited);

dist_tagged = histc(activation_peak(tagged),edges);
dist_tagged = dist_tagged(1:end-1);
dist_tagged = dist_tagged / sum(dist_tagged);

figure
plot(centers,dist_activated,'Color','r','LineWidth',3)
hold on
plot(centers,dist_inhibited,'Color','b','LineWidth',3)
plot(centers,dist_tagged,'Color','k','LineWidth',3)

%% cdf #3

edges = -0.5:1:100.5;
centers = (edges(1:end-1) + edges(2:end)) / 2;
dist_activated = histc(activation_peak(activated),edges);
dist_activated = dist_activated(1:end-1);
dist_activated = dist_activated / sum(dist_activated);

dist_inhibited = histc(inhibition_peak(inhibited),edges);
dist_inhibited = dist_inhibited(1:end-1);
dist_inhibited = dist_inhibited / sum(dist_inhibited);

dist_tagged = histc(activation_peak(tagged),edges);
dist_tagged = dist_tagged(1:end-1);
dist_tagged = dist_tagged / sum(dist_tagged);

figure
stairs(centers,cumsum(dist_activated),'Color',[0.32 0.19 0.19],'LineWidth',3)
hold on
stairs(centers,cumsum(dist_inhibited),'Color',[0.48 0.06 0.89],'LineWidth',3)
stairs(centers,cumsum(dist_tagged),'Color',[0 0.7 0],'LineWidth',3)

set(gca,'box','off','XTick',0:10:70,'YTick',[0:0.5:1],...
    'FontSize',16,'TickDir','out','XLim',[0 71],'Ylim',[0 1.05])
xlabel('Peak time')