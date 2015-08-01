function vipplotresponse2_A1_fcn
%VIPPLOTRESPONSE2_A1_FCN   Peak effect vs. firing rate change.
%   VIPPLOTRESPONSE2_A1_FCN plots inhibition/activation peak time from
%   light-stimulus aligned PSTHs vs relative firing rate change (max/min to
%   baseline) on log scale. Only significantly influenced cells are
%   included (p<0.05, two sided MW U-test; see VIPISINFLUENCED3B). Range of
%   half-peak crossings is indicated by horizontal whiskers. Cells
%   inhibited to zero firing rate are plotted below all other points, with
%   y value randomly jittered to prevent overlap. Cells are separated to
%   directly activated ('tagged'), inhibited and delayed activated
%   (referred to as 'activated') groups. The peak time pdf for the
%   different groups is superimposed. The pdf for the activated group is
%   smoothed by using overlapping bins. The relative scaling of pdfs is
%   arbitrary. Peak time cdf is also plotted.
%
%   The function also plots spike width vs. firing rate, spike width cdf
%   and firing rate cdf. For these cdfs, the 'all cells' group (gray)
%   contains only those cells that are above the lower limit for possible
%   inhibition detection (1 Hz).
%
%   See also VIPPLOTRESPONSE2_A1 and VIPISINFLUENCED_3B.

% Load PSTH variables
% load('C:\Balazs\_analysis\VIP\A1_psth_variables_revision1.mat')
load('C:\Balazs\_analysis\VIP\A1_psth_variables.mat')

% Groups
frlim = 1;   % lower firing rate limit for detecting inhibition
tagged = logical(vldty) & (isact==2);   % tagged cells
inx_act = logical(vldty) & (isact==1) & (isact~=2);  % activated cells
inx_inh = logical(vldty) & (isinh) & (baseline>frlim) & (isact~=2);   % inhibited cells; firing rate > 1Hz criterion
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

% Relative FR change of inhibited cells
finh = minvalue(inhibited) ./ baseline(inhibited);   % rel. FR change for inhibited cells
r = rand(1,sum(finh==0)) / 700;
finh(finh==0) = r + 0.001;   % jitter cells that are inhibited to FR=0

% Define colors
cact = [1 0.55 0.33];   % activated cells
cinh = [0.48 0.06 0.89];   % inhibited cells
cina = [1 0 1];   % inhibited-activated cells
ctag = [0 0.7 0];   % tagged cells
cact_line = [1 0.8 0.4];   % whiskers for activated cells
cinh_line = [0.85 0.7 1];   % whiskers for inhibited cells
ctag_line = [0.85 0.7 1];   % whiskers for tagged cells

% Plot peak time vs relative FR change
figure
semilogy(activation_peak(activated),maxvalue(activated)./baseline(activated),'o','MarkerEdgeColor',cact,'MarkerFaceColor',cact,'MarkerSize',8)
hold on   % log plot should precede the lines, otherwise the plot is screwed
semilogy(inhibition_peak(inhibited),finh,'o','MarkerEdgeColor',cinh,'MarkerFaceColor',cinh,'MarkerSize',8)
semilogy(activation_peak(tagged),maxvalue(tagged)./baseline(tagged),'o','MarkerEdgeColor',ctag,'MarkerFaceColor',ctag,'MarkerSize',8)
semilogy(inhibition_peak(inhibited(iainx)),finh(iainx),'o','MarkerEdgeColor',cina,'MarkerFaceColor',cina,'MarkerSize',8)

% Plot whiskers for half-peak ranges
% for k = 1:length(activated)   % activated
%     line([activation_start(activated(k)) activation_end(activated(k))],...
%         [maxvalue(activated(k))./baseline(activated(k)) maxvalue(activated(k))./baseline(activated(k))],'Color',cact_line)
% end
% for k = 1:length(tagged)   % tagged
%     line([activation_start(tagged(k)) activation_end(tagged(k))],...
%         [maxvalue(tagged(k))./baseline(tagged(k)) maxvalue(tagged(k))./baseline(tagged(k))],'Color',ctag_line)
% end
% for k = 1:length(inhibited)   % inhibited
%     line([inhibition_start(inhibited(k)) inhibition_end(inhibited(k))],...
%         [finh(k) finh(k)],'Color',cinh_line)
% end
% for k = 1:length(inhibited_activated)
%     line([inhibition_start(inhibited_activated(k)) inhibition_end(inhibited_activated(k))],...
%         [finh(k) finh(k)],'Color',cinh_line)
% end

% Scatter plot
semilogy(activation_peak(activated),maxvalue(activated)./baseline(activated),'o','MarkerEdgeColor',cact,'MarkerFaceColor',cact,'MarkerSize',8)
hold on   % plot dots again to overlay them on the lines
semilogy(inhibition_peak(inhibited),finh,'o','MarkerEdgeColor',cinh,'MarkerFaceColor',cinh,'MarkerSize',8)
semilogy(activation_peak(tagged),maxvalue(tagged)./baseline(tagged),'o','MarkerEdgeColor',ctag,'MarkerFaceColor',ctag,'MarkerSize',8)
semilogy(inhibition_peak(inhibited(iainx)),finh(iainx),'o','MarkerEdgeColor',cina,'MarkerFaceColor',cina,'MarkerSize',8)

% Format figure
line([0 200],[1 1],'Color','k')
set(gca,'box','off','XTick',0:50:200,'YTick',[0.01 0.1 1 10 100 1000],'YMinorTick','off',...
    'FontSize',16,'TickDir','out','XLim',[0 201])
xlabel('Peak time')
ylabel('Relative change in firing rate')

% Peak time pdf for activated cells
edges = [0:5:170; 30:5:200];     % histogram bin limits
centers_activated = (edges(1,:) + edges(2,:)) / 2;     % histogram bin centers
dist_activated = histc_overlap(activation_peak(activated),edges);   % histogram with overlapping bis
dist_activated = dist_activated ./ (edges(2,:) - edges(1,:));   % normalize to bin size
dist_activated = dist_activated / sum(dist_activated);   % normalize sum to 1

% Peak time pdf for tagged cells
edges = [0 0.7 1.5 2.5 7 10 20 30 60 90 120 150 180 200];     % histogram bin limits
centers_tagged = (edges(1:end-1) + edges(2:end)) / 2;     % histogram bin centers
dist_tagged = histc(activation_peak(tagged),edges);   % histogram with non-overlapping bis
dist_tagged = dist_tagged(1:end-1) ./ diff(edges');   % normalize to bin size
dist_tagged = dist_tagged / sum(dist_tagged);   % normalize sum to 1

% Peak time pdf for inhibited cells
edges = [0 2 7 11 16 20 30 60 90 120 150 180 200];     % histogram bin limits
centers_inhibited = (edges(1:end-1) + edges(2:end)) / 2;     % histogram bin centers
dist_inhibited = histc(inhibition_peak(inhibited),edges);   % histogram with non-overlapping bis
dist_inhibited = dist_inhibited(1:end-1) ./ diff(edges');   % normalize to bin size
dist_inhibited = dist_inhibited / sum(dist_inhibited);   % normalize sum to 1

% Overlay pdfs
dist_activated2 = exp(dist_activated * 60);  % rescale activated pdf and convert it for the log plot
dist_inhibited2 = exp(-dist_inhibited * 17);  % rescale inhibited pdf and convert it for the log plot
dist_tagged2 = exp(dist_tagged * 8);  % rescale tagged pdf and convert it for the log plot
plot([0 centers_activated],[1 dist_activated2],'Color',cact,'LineWidth',3)   % activated pdf
hold on
plot(centers_inhibited,dist_inhibited2,'Color',cinh,'LineWidth',3)   % inhibited pdf
plot(centers_tagged,dist_tagged2,'Color',ctag,'LineWidth',3)   % tagged pdf

% keyboard

% Peak time cdf
edges = -0.5:1:200.5;   % histogram bin limits
centers = (edges(1:end-1) + edges(2:end)) / 2;   % histogram bin centers
dist_activated = histc(activation_peak(activated),edges);   % activated cells
dist_activated = dist_activated(1:end-1);
dist_activated = dist_activated / sum(dist_activated);

dist_inhibited = histc(inhibition_peak(inhibited),edges);   % inhibited cells
dist_inhibited = dist_inhibited(1:end-1);
dist_inhibited = dist_inhibited / sum(dist_inhibited);

dist_tagged = histc(activation_peak(tagged),edges);   % tagged cells
dist_tagged = dist_tagged(1:end-1);
dist_tagged = dist_tagged / sum(dist_tagged);

% Plot peak time cdf
figure
stairs(centers,cumsum(dist_tagged),'Color',ctag,'LineWidth',3)
hold on
stairs(centers,cumsum(dist_inhibited),'Color',cinh,'LineWidth',3)
stairs(centers,cumsum(dist_activated),'Color',cact,'LineWidth',3)

% Format figure
set(gca,'box','off','XTick',0:50:200,'YTick',[0:0.5:1],...
    'FontSize',16,'TickDir','out','XLim',[0 201],'Ylim',[0 1.05])
xlabel('Peak time')

% Plot Spike width vs. firing rate
g2inx = baseline > frlim;   % indices for cells above firing rate limit for inhibtion
s2inx = baseline <= frlim;   % indices for cells below firing rate limit for inhibtion
figure;
plot(spike_width(g2inx),baseline(g2inx),'o','MarkerEdgeColor',[0.31 0.31 0.31],'MarkerSize',3,'MarkerFaceColor',[0.31 0.31 0.31])  % unidentified cells >1Hz
hold on
plot(spike_width(s2inx),baseline(s2inx),'o','Color',[0.7 0.7 0.7],'MarkerSize',3,'MarkerFaceColor',[0.7 0.7 0.7])  % unidentified cells <=1Hz
plot(spike_width(activated),baseline(activated),'o','MarkerEdgeColor',cact,'MarkerFaceColor',cact,'MarkerSize',8)  % activated cells
plot(spike_width(inhibited),baseline(inhibited),'o','MarkerEdgeColor',cinh,'MarkerFaceColor',cinh,'MarkerSize',8)  % inhibited cells
plot(spike_width(inhibited_activated),baseline(inhibited_activated),'o','MarkerEdgeColor',cina,'MarkerFaceColor',cina,'MarkerSize',8)  % inhibited-activated cells
plot(spike_width(tagged),baseline(tagged),'o','MarkerEdgeColor',ctag,'MarkerFaceColor',ctag,'MarkerSize',8)  % tagged cells

% Format figure
line([100 600],[frlim frlim],'Color','k','LineStyle',':')
% plot(spike_width(both),baseline(both),'go','MarkerfaceColor','g')
set(gca,'box','off','XTick',0:200:600,'YTick',0:10:50,'FontSize',16,'TickDir','out')
xlabel('Spike width')
ylabel('Firing rate')

% Spike width cdf
edges = -0.5:1:600.5;   % histogram bin limits
centers = (edges(1:end-1) + edges(2:end)) / 2;   % histogram bin centers
dist_activated = histc(spike_width(activated),edges);
dist_activated = histc(spike_width(intersect(activated,find(baseline>frlim))),edges);   % activated cells
dist_activated = dist_activated(1:end-1);
dist_activated = dist_activated / sum(dist_activated);

dist_inhibited = histc(spike_width(inhibited),edges);   % inhibited cells
dist_inhibited = dist_inhibited(1:end-1);
dist_inhibited = dist_inhibited / sum(dist_inhibited);

dist_all = histc(spike_width(baseline>frlim),edges);   % all cells above firing rate limit  
dist_all = dist_all(1:end-1);
dist_all = dist_all / sum(dist_all);

% Plot spike width cdf
figure
stairs(centers,cumsum(dist_activated),'Color',cact,'LineWidth',3)   % activated cells
hold on
stairs(centers,cumsum(dist_inhibited),'Color',[0.48 0.06 0.89],'LineWidth',3)   % inhibited cells
stairs(centers,cumsum(dist_all),'Color',[0.31 0.31 0.31],'LineWidth',3)   % all cells above firing rate limit

% Format figure
set(gca,'box','off','XTick',0:200:600,'YTick',0:0.5:1,'FontSize',16,'TickDir','out','YLim',[0 1.05])
xlabel('Spike width')

% Firing rate cdf
edges = 0:0.1:50;   % histogram bin limits
centers = (edges(1:end-1) + edges(2:end)) / 2;   % histogram bin edges
dist_activated = histc(baseline(activated),edges);
dist_activated = histc(baseline(intersect(activated,find(baseline>frlim))),edges);   % activated cells
dist_activated = dist_activated(1:end-1);
dist_activated = dist_activated / sum(dist_activated);

dist_inhibited = histc(baseline(inhibited),edges);   % inhibited cells
dist_inhibited = dist_inhibited(1:end-1);
dist_inhibited = dist_inhibited / sum(dist_inhibited);

dist_all = histc(baseline(baseline>frlim),edges);   % all cells above firing rate limit
dist_all = dist_all(1:end-1);
dist_all = dist_all / sum(dist_all);

% Plot firing rate cdf
figure
plot(centers,cumsum(dist_activated),'Color',cact,'LineWidth',3)   % activated cells
hold on
plot(centers,cumsum(dist_inhibited),'Color',cinh,'LineWidth',3)   % inhibited cells
plot(centers,cumsum(dist_all),'Color',[0.31 0.31 0.31],'LineWidth',3)   % all cells above firing rate limit 

% Format figure
set(gca,'box','off','XTick',0:10:30,'YTick',0:0.5:1,'FontSize',16,'TickDir','out','YLim',[0 1.05],'Xlim',[0 30])
xlabel('Firing rate')

% -------------------------------------------------------------------------
function mn = histc_overlap(x,edges)
%HISTC_OVERLAP   Histogram with overlapping bins.

% Convert edges to 2xN
[ne me] = size(edges);
if ne > me
    edges = edges';
    [~, me] = size(edges);
end

% Calculate counts
n = me;   % number of bins
mn = nan(1,n);
for k = 1:n
    mn(k) = sum((x>edges(1,k)&x<edges(2,k)));  % count
end