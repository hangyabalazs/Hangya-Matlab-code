%% load

load('C:\Balazs\_analysis\VIP\psth_variables.mat')

%% groups

tagged = logical(vldty) & (isact==2);
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
inhibited_short = intersect(inhibited,find(inhibition_peak<40));
inhibited_long = intersect(inhibited,find(inhibition_peak>=40));

%% Plots

figure;
g2inx = baseline > 2;
s2inx = baseline <= 2;
plot(spike_width(g2inx),baseline(g2inx),'o','MarkerEdgeColor',[0.31 0.31 0.31],'MarkerSize',3,'MarkerFaceColor',[0.31 0.31 0.31])
hold on
plot(spike_width(s2inx),baseline(s2inx),'o','Color',[0.7 0.7 0.7],'MarkerSize',3,'MarkerFaceColor',[0.7 0.7 0.7])
plot(spike_width(activated),baseline(activated),'o','MarkerEdgeColor',[0.32 0.19 0.19],'MarkerFaceColor',[0.32 0.19 0.19],'MarkerSize',8)
plot(spike_width(inhibited_short),baseline(inhibited_short),'o','MarkerEdgeColor',[0.48 0.06 0.89],'MarkerFaceColor',[0.48 0.06 0.89],'MarkerSize',8)
plot(spike_width(tagged),baseline(tagged),'o','MarkerEdgeColor',[0 0.7 0],'MarkerFaceColor',[0 0.7 0],'MarkerSize',8)
line([100 600],[2 2],'Color','k','LineStyle',':')
% plot(spike_width(both),baseline(both),'go','MarkerfaceColor','g')
set(gca,'box','off','XTick',100:100:500,'YTick',0:10:50,'FontSize',16,'TickDir','out')
xlabel('Spike width')
ylabel('Firing rate')

%% cdf - spike width

edges = -0.5:1:600.5;
centers = (edges(1:end-1) + edges(2:end)) / 2;
dist_activated = histc(spike_width(activated),edges);
dist_activated = dist_activated(1:end-1);
dist_activated = dist_activated / sum(dist_activated);

dist_inhibited = histc(spike_width(inhibited_short),edges);
dist_inhibited = dist_inhibited(1:end-1);
dist_inhibited = dist_inhibited / sum(dist_inhibited);

dist_all = histc(spike_width(baseline>2),edges);
dist_all = dist_all(1:end-1);
dist_all = dist_all / sum(dist_all);

figure
stairs(centers,cumsum(dist_activated),'Color',[0.32 0.19 0.19],'LineWidth',3)
hold on
stairs(centers,cumsum(dist_inhibited),'Color',[0.48 0.06 0.89],'LineWidth',3)
stairs(centers,cumsum(dist_all),'Color',[0.31 0.31 0.31],'LineWidth',3)
set(gca,'box','off','XTick',100:100:500,'FontSize',16,'TickDir','out','YLim',[0 1.05])
xlabel('Spike width')

%% colors

cact = [1 0.55 0.33];   % activated cells
cinh = [0.48 0.06 0.89];   % inhibited cells
cina = [1 0 1];   % inhibited-activated cells
ctag = [0 0.7 0];   % tagged cells

%% cdf - FR

edges = 0:0.1:50;
centers = (edges(1:end-1) + edges(2:end)) / 2;
dist_activated = histc(baseline(activated),edges);
dist_activated = dist_activated(1:end-1);
dist_activated = dist_activated / sum(dist_activated);

dist_inhibited = histc(baseline(inhibited_short),edges);
dist_inhibited = dist_inhibited(1:end-1);
dist_inhibited = dist_inhibited / sum(dist_inhibited);

dist_all = histc(baseline(baseline>2),edges);
dist_all = dist_all(1:end-1);
dist_all = dist_all / sum(dist_all);

figure
plot(centers,cumsum(dist_activated),'Color',cact,'LineWidth',3)
hold on
plot(centers,cumsum(dist_inhibited),'Color',cinh,'LineWidth',3)
plot(centers,cumsum(dist_all),'Color',[0.31 0.31 0.31],'LineWidth',3)