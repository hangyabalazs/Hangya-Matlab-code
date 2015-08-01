function vipplotresponse2_A1_fcn

% load

load('C:\Balazs\_analysis\VIP\A1_psth_variables.mat')

% groups

tagged = logical(vldty) & (isact==2);
% inx_act = logical(vldty) & (isact==1) & (baseline>1);
inx_act = logical(vldty) & (isact==1);
inx_inh = logical(vldty) & (isinh) & (baseline>1) & (isact~=2);
activated = setdiff(find(inx_act&~inx_inh),tagged);
inhibited = find(inx_inh&~inx_act);
ai = find(inx_act&inx_inh);
inx = activation_peak(ai) < inhibition_peak(ai);
activated_inhibited = sort(ai(inx));
inhibited_activated = ai(~inx);
activated = [activated; activated_inhibited];
inhibited = [inhibited; inhibited_activated];
tagged = find(tagged);

% peak time - with lines, -Inf added, - with pdf

finh = minvalue(inhibited) ./ baseline(inhibited);
r = rand(1,sum(finh==0)) / 500;
finh(finh==0) = r + 0.002;

figure
semilogy(activation_peak(activated),maxvalue(activated)./baseline(activated),'o','MarkerEdgeColor',[1 0.55 0.33],'MarkerFaceColor',[1 0.55 0.33],'MarkerSize',8)
hold on
semilogy(inhibition_peak(inhibited),finh,'o','MarkerEdgeColor',[0.48 0.06 0.89],'MarkerFaceColor',[0.48 0.06 0.89],'MarkerSize',8)
semilogy(activation_peak(tagged),maxvalue(tagged)./baseline(tagged),'o','MarkerEdgeColor',[0 0.7 0],'MarkerFaceColor',[0 0.7 0],'MarkerSize',8)

% for k = 1:length(activated)
%     line([activation_start(activated(k)) activation_end(activated(k))],...
%         [maxvalue(activated(k))./baseline(activated(k)) maxvalue(activated(k))./baseline(activated(k))],'Color',[1 0.8 0.4])
% end
% for k = 1:length(tagged)
%     line([activation_start(tagged(k)) activation_end(tagged(k))],...
%         [maxvalue(tagged(k))./baseline(tagged(k)) maxvalue(tagged(k))./baseline(tagged(k))],'Color',[0 0.7 0])
% end
% for k = 1:length(inhibited)
%     line([inhibition_start(inhibited(k)) inhibition_end(inhibited(k))],...
%         [finh(k) finh(k)],'Color',[0.85 0.7 1])
% end

semilogy(activation_peak(activated),maxvalue(activated)./baseline(activated),'o','MarkerEdgeColor',[1 0.55 0.33],'MarkerFaceColor',[1 0.55 0.33],'MarkerSize',8)
semilogy(inhibition_peak(inhibited),finh,'o','MarkerEdgeColor',[0.48 0.06 0.89],'MarkerFaceColor',[0.48 0.06 0.89],'MarkerSize',8)
semilogy(activation_peak(tagged),maxvalue(tagged)./baseline(tagged),'o','MarkerEdgeColor',[0 0.7 0],'MarkerFaceColor',[0 0.7 0],'MarkerSize',8)
line([0 200],[1 1],'Color','k')
% set(gca,'box','off','XTick',0:50:200,'YTick',[0.01 0.1 1 10 100 1000],'YMinorTick','off',...
%     'FontSize',16,'TickDir','out','XLim',[0 201],'Ylim',[0.001 1000])
xlabel('Peak time')
ylabel('Relative change in firing rate')

% edges = [0 15:20:150 180 200];
% centers_activated = (edges(1:end-1) + edges(2:end)) / 2;
edges = [0:5:170; 30:5:200];     % histogram bin limits
centers_activated = (edges(1,:) + edges(2,:)) / 2;     % histogram bin centers
dist_activated = histc_overlap(activation_peak(activated),edges);
% dist_activated = dist_activated(1:end-1) ./ diff(edges');
dist_activated = dist_activated ./ (edges(2,:) - edges(1,:));
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
dist_activated2 = exp(dist_activated * 60);
dist_inhibited2 = exp(-dist_inhibited * 15);
dist_tagged2 = exp(dist_tagged * 10);
% plot(centers_activated,dist_activated2,'Color',[1 0.55 0.33],'LineWidth',3)
plot([0 centers_activated],[1 dist_activated2],'Color',[1 0.55 0.33],'LineWidth',3)
hold on
plot(centers_inhibited,dist_inhibited2,'Color',[0.48 0.06 0.89],'LineWidth',3)
plot(centers_tagged,dist_tagged2,'Color',[0 0.7 0],'LineWidth',3)

keyboard


% -------------------------------------------------------------------------
function mn = histc_overlap(x,edges)

[ne me] = size(edges);
if ne > me
    edges = edges';
    [ne me] = size(edges);
end

% Calculate conditional mean values
n = me;
mn = zeros(1,n);
for k = 1:n
    mn(k) = sum((x>edges(1,k)&x<edges(2,k)));
end