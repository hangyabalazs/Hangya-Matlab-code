%% load

load('C:\Balazs\_analysis\VIP\A1_psth_variables.mat')

%% groups

tagged = logical(vldty) & (isact==2);   % tagged cells
inx_act = logical(vldty) & (isact==1) & (isact~=2);  % activated cells
inx_inh = logical(vldty) & (isinh) & (baseline>1) & (isact~=2);   % inhibited cells; firing rate > 1Hz criterion
activated = find(inx_act&~inx_inh);  % indices of activated only cells
inhibited = find(inx_inh&~inx_act);  % indices of inhibited only cells
ai = find(inx_act&inx_inh);   % indices of cells with dual effect
inx = activation_peak(ai) < inhibition_peak(ai);
activated_inhibited = sort(ai(inx));  % activated, then inhibited
inhibited_activated = ai(~inx);   % inhibited, then activated
activated = [activated; activated_inhibited];   % activated and activated-inhibited
inhibited = [inhibited; inhibited_activated];   % inhibited and inhibited-activated
tagged = find(tagged);   % indices of tagged cells
[~, iainx] = intersect(inhibited,inhibited_activated);   % indices of inh-act within inh

%% peak time

figure
semilogy(activation_peak(activated),maxvalue(activated)./baseline(activated),'o','MarkerEdgeColor',[1 0.55 0.33],'MarkerFaceColor',[1 0.55 0.33],'MarkerSize',8)
hold on
semilogy(inhibition_peak(inhibited),minvalue(inhibited)./baseline(inhibited),'o','MarkerEdgeColor',[0.48 0.06 0.89],'MarkerFaceColor',[0.48 0.06 0.89],'MarkerSize',8)
semilogy(activation_peak(tagged),maxvalue(tagged)./baseline(tagged),'o','MarkerEdgeColor',[0 0.7 0],'MarkerFaceColor',[0 0.7 0],'MarkerSize',8)
line([0 200],[1 1],'Color','k')
set(gca,'box','off','XTick',0:20:200,'YTick',[0.001 0.01 0.1 1 10 100 1000],'YMinorTick','off',...
    'FontSize',16,'TickDir','out','XLim',[0 201],'Ylim',[0.01 1000])
xlabel('Peak time')
ylabel('Relative change in firing rate')

%% peak time - with lines

figure
semilogy(activation_peak(activated),maxvalue(activated)./baseline(activated),'o','MarkerEdgeColor',[1 0.55 0.33],'MarkerFaceColor',[1 0.55 0.33],'MarkerSize',8)
hold on
semilogy(inhibition_peak(inhibited),minvalue(inhibited)./baseline(inhibited),'o','MarkerEdgeColor',[0.48 0.06 0.89],'MarkerFaceColor',[0.48 0.06 0.89],'MarkerSize',8)
semilogy(activation_peak(tagged),maxvalue(tagged)./baseline(tagged),'o','MarkerEdgeColor',[0 0.7 0],'MarkerFaceColor',[0 0.7 0],'MarkerSize',8)
for k = 1:length(activated)
    line([activation_start(activated(k)) activation_end(activated(k))],...
        [maxvalue(activated(k))./baseline(activated(k)) maxvalue(activated(k))./baseline(activated(k))],'Color',[1 0.55 0.33])
end
for k = 1:length(tagged)
    line([activation_start(tagged(k)) activation_end(tagged(k))],...
        [maxvalue(tagged(k))./baseline(tagged(k)) maxvalue(tagged(k))./baseline(tagged(k))],'Color',[0 0.7 0])
end
line([0 200],[1 1],'Color','k')
set(gca,'box','off','XTick',0:50:200,'YTick',[0.001 0.01 0.1 1 10 100 1000],'YMinorTick','off',...
    'FontSize',16,'TickDir','out','XLim',[0 201],'Ylim',[0.01 1000])
xlabel('Peak time')
ylabel('Relative change in firing rate')

%% peak time - with lines, -Inf added

finh = minvalue(inhibited) ./ baseline(inhibited);
r = rand(1,sum(finh==0)) / 500;
finh(finh==0) = r + 0.002;

figure
semilogy(activation_peak(activated),maxvalue(activated)./baseline(activated),'o','MarkerEdgeColor',[1 0.55 0.33],'MarkerFaceColor',[1 0.55 0.33],'MarkerSize',8)
hold on
semilogy(inhibition_peak(inhibited),finh,'o','MarkerEdgeColor',[0.48 0.06 0.89],'MarkerFaceColor',[0.48 0.06 0.89],'MarkerSize',8)
semilogy(activation_peak(tagged),maxvalue(tagged)./baseline(tagged),'o','MarkerEdgeColor',[0 0.7 0],'MarkerFaceColor',[0 0.7 0],'MarkerSize',8)

for k = 1:length(activated)
    line([activation_start(activated(k)) activation_end(activated(k))],...
        [maxvalue(activated(k))./baseline(activated(k)) maxvalue(activated(k))./baseline(activated(k))],'Color',[1 0.8 0.4])
end
for k = 1:length(tagged)
    line([activation_start(tagged(k)) activation_end(tagged(k))],...
        [maxvalue(tagged(k))./baseline(tagged(k)) maxvalue(tagged(k))./baseline(tagged(k))],'Color',[0 0.7 0])
end
for k = 1:length(inhibited)
    line([inhibition_start(inhibited(k)) inhibition_end(inhibited(k))],...
        [finh(k) finh(k)],'Color',[0.85 0.7 1])
end

semilogy(activation_peak(activated),maxvalue(activated)./baseline(activated),'o','MarkerEdgeColor',[1 0.55 0.33],'MarkerFaceColor',[1 0.55 0.33],'MarkerSize',8)
semilogy(inhibition_peak(inhibited),finh,'o','MarkerEdgeColor',[0.48 0.06 0.89],'MarkerFaceColor',[0.48 0.06 0.89],'MarkerSize',8)
semilogy(activation_peak(tagged),maxvalue(tagged)./baseline(tagged),'o','MarkerEdgeColor',[0 0.7 0],'MarkerFaceColor',[0 0.7 0],'MarkerSize',8)
line([0 200],[1 1],'Color','k')
set(gca,'box','off','XTick',0:50:200,'YTick',[0.01 0.1 1 10 100 1000],'YMinorTick','off',...
    'FontSize',16,'TickDir','out','XLim',[0 201],'Ylim',[0.001 1000])
xlabel('Peak time')
ylabel('Relative change in firing rate')

%% peak time - with lines, -Inf added, - with pdf

finh = minvalue(inhibited) ./ baseline(inhibited);
r = rand(1,sum(finh==0)) / 500;
finh(finh==0) = r + 0.002;

figure
semilogy(activation_peak(activated),maxvalue(activated)./baseline(activated),'o','MarkerEdgeColor',[1 0.55 0.33],'MarkerFaceColor',[1 0.55 0.33],'MarkerSize',8)
hold on
semilogy(inhibition_peak(inhibited),finh,'o','MarkerEdgeColor',[0.48 0.06 0.89],'MarkerFaceColor',[0.48 0.06 0.89],'MarkerSize',8)
semilogy(activation_peak(tagged),maxvalue(tagged)./baseline(tagged),'o','MarkerEdgeColor',[0 0.7 0],'MarkerFaceColor',[0 0.7 0],'MarkerSize',8)

for k = 1:length(activated)
    line([activation_start(activated(k)) activation_end(activated(k))],...
        [maxvalue(activated(k))./baseline(activated(k)) maxvalue(activated(k))./baseline(activated(k))],'Color',[1 0.8 0.4])
end
for k = 1:length(tagged)
    line([activation_start(tagged(k)) activation_end(tagged(k))],...
        [maxvalue(tagged(k))./baseline(tagged(k)) maxvalue(tagged(k))./baseline(tagged(k))],'Color',[0 0.7 0])
end
for k = 1:length(inhibited)
    line([inhibition_start(inhibited(k)) inhibition_end(inhibited(k))],...
        [finh(k) finh(k)],'Color',[0.85 0.7 1])
end

semilogy(activation_peak(activated),maxvalue(activated)./baseline(activated),'o','MarkerEdgeColor',[1 0.55 0.33],'MarkerFaceColor',[1 0.55 0.33],'MarkerSize',8)
semilogy(inhibition_peak(inhibited),finh,'o','MarkerEdgeColor',[0.48 0.06 0.89],'MarkerFaceColor',[0.48 0.06 0.89],'MarkerSize',8)
semilogy(activation_peak(tagged),maxvalue(tagged)./baseline(tagged),'o','MarkerEdgeColor',[0 0.7 0],'MarkerFaceColor',[0 0.7 0],'MarkerSize',8)
line([0 200],[1 1],'Color','k')
set(gca,'box','off','XTick',0:50:200,'YTick',[0.01 0.1 1 10 100 1000],'YMinorTick','off',...
    'FontSize',16,'TickDir','out','XLim',[0 201],'Ylim',[0.001 1000])
xlabel('Peak time')
ylabel('Relative change in firing rate')

edges = [0 15:20:150 180 200];
centers_activated = (edges(1:end-1) + edges(2:end)) / 2;
dist_activated = histc(activation_peak(activated),edges);
dist_activated = dist_activated(1:end-1) ./ diff(edges');
dist_activated = dist_activated / sum(dist_activated);

edges = [0 0.7 1.5 2.5 4 7 10 20 30 60 90 120 150 180 200];
centers_tagged = (edges(1:end-1) + edges(2:end)) / 2;
dist_tagged = histc(activation_peak(tagged),edges);
dist_tagged = dist_tagged(1:end-1) ./ diff(edges');
dist_tagged = dist_tagged / sum(dist_tagged);

edges = [0 2 5 10 15 20 30 60 90 120 150 180 200];
centers_inhibited = (edges(1:end-1) + edges(2:end)) / 2;
dist_inhibited = histc(inhibition_peak(inhibited),edges);
dist_inhibited = dist_inhibited(1:end-1) ./ diff(edges');
dist_inhibited = dist_inhibited / sum(dist_inhibited);

% figure
dist_activated2 = exp(dist_activated * 15);
dist_inhibited2 = exp(-dist_inhibited * 15);
dist_tagged2 = exp(dist_tagged * 10);
plot(centers_activated,dist_activated2,'Color',[1 0.55 0.33],'LineWidth',3)
hold on
plot(centers_inhibited,dist_inhibited2,'Color',[0.48 0.06 0.89],'LineWidth',3)
plot(centers_tagged,dist_tagged2,'Color',[0 0.7 0],'LineWidth',3)


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

edges = -0.5:1:200.5;
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
stairs(centers,cumsum(dist_tagged),'Color',[0 0.7 0],'LineWidth',3)
hold on
stairs(centers,cumsum(dist_inhibited),'Color',[0.48 0.06 0.89],'LineWidth',3)
stairs(centers,cumsum(dist_activated),'Color',[1 0.55 0.33],'LineWidth',3)

set(gca,'box','off','XTick',0:50:200,'YTick',[0:0.5:1],...
    'FontSize',16,'TickDir','out','XLim',[0 201],'Ylim',[0 1.05])
xlabel('Peak time')

%% sound responsive and non-responsive activated cells

% Groups
Sound = [46 56 162 167 168 169 99 126 28 146];
NoSound = [47 49 57 60 63 65 66 149 177 173];

figure
plot(maxvalue(Sound)./baseline(Sound),activation_time(Sound),'ro','MarkerFaceColor','r')
hold on
plot(maxvalue(NoSound)./baseline(NoSound),activation_time(NoSound),'o','MarkerFaceColor','b')