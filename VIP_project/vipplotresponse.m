%% load

load('C:\Balazs\_analysis\VIP\vipisinfluenced11\validity.mat')
load('C:\Balazs\_analysis\VIP\vipisinfluenced11\p_val.mat')

%% without FR criterion

inx_act = logical(vldty') & (p_act<0.015);
inx_inh = logical(vldty') & (p_inh<0.015);

%% with FR criterion

inx_act = logical(vldty') & (p_act<0.015) & (baseline>2);
inx_act(105) = 1;   % tagged
inx_inh = logical(vldty') & (p_inh<0.015) & (baseline>2);

%% peak time

figure
plot(activation_peak(inx_act),maxvalue(inx_act)-baseline(inx_act),'o','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',10)
hold on
plot(inhibition_peak(inx_inh),minvalue(inx_inh)-baseline(inx_inh),'o','MarkerFaceColor','b','MarkerEdgeColor','b','MarkerSize',10)


figure
plot(activation_peak(inx_act),maxvalue(inx_act)./baseline(inx_act),'o','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',10)
hold on
plot(inhibition_peak(inx_inh),minvalue(inx_inh)./baseline(inx_inh),'o','MarkerFaceColor','b','MarkerEdgeColor','b','MarkerSize',10)


figure
plot(activation_peak(inx_act),log10(maxvalue(inx_act)./baseline(inx_act)),'o','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',10)
hold on
plot(inhibition_peak(inx_inh),log10(minvalue(inx_inh)./baseline(inx_inh)),'o','MarkerFaceColor','b','MarkerEdgeColor','b','MarkerSize',10)

%% colors

figure
plot(activation_peak(inx_act),log10(maxvalue(inx_act)./baseline(inx_act)),'o','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',10)
hold on
plot(inhibition_peak(inx_inh),log10(minvalue(inx_inh)./baseline(inx_inh)),'o','MarkerFaceColor','b','MarkerEdgeColor','b','MarkerSize',10)
plot(inhibition_peak(inx_act&inx_inh),log10(minvalue(inx_act&inx_inh)./baseline(inx_act&inx_inh)),'o','MarkerFaceColor','m','MarkerEdgeColor','m','MarkerSize',10)
plot(activation_peak(inx_act&inx_inh),log10(maxvalue(inx_act&inx_inh)./baseline(inx_act&inx_inh)),'o','MarkerFaceColor','y','MarkerEdgeColor','y','MarkerSize',10)

%% start time

figure
plot(activation_start(inx_act),maxvalue(inx_act)-baseline(inx_act),'o','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',10)
hold on
plot(inhibition_start(inx_inh),minvalue(inx_inh)-baseline(inx_inh),'o','MarkerFaceColor','b','MarkerEdgeColor','b','MarkerSize',10)


figure
plot(activation_start(inx_act),maxvalue(inx_act)./baseline(inx_act),'o','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',10)
hold on
plot(inhibition_start(inx_inh),minvalue(inx_inh)./baseline(inx_inh),'o','MarkerFaceColor','b','MarkerEdgeColor','b','MarkerSize',10)


figure
plot(activation_start(inx_act),log10(maxvalue(inx_act)./baseline(inx_act)),'o','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',10)
hold on
plot(inhibition_start(inx_inh),log10(minvalue(inx_inh)./baseline(inx_inh)),'o','MarkerFaceColor','b','MarkerEdgeColor','b','MarkerSize',10)
plot(inhibition_start(inx_act&inx_inh),log10(minvalue(inx_act&inx_inh)./baseline(inx_act&inx_inh)),'o','MarkerFaceColor','m','MarkerEdgeColor','m','MarkerSize',10)
plot(activation_start(inx_act&inx_inh),log10(maxvalue(inx_act&inx_inh)./baseline(inx_act&inx_inh)),'o','MarkerFaceColor','y','MarkerEdgeColor','y','MarkerSize',10)


%% pdf

tagged = find(inx_act&activation_start<3);
activated = setdiff(find(inx_act),tagged);
inhibited = find(inx_inh);

edges = 0:5:100;
edges = -2.5:5:102.5;
% edges = 0:3:100;
% edges = [0 1 2 4 8 16 32 64 128];
% edges = 0:10:100;
% edges = [0 3.^(0:5)];
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

%% cdf

figure
plot(centers,cumsum(dist_activated),'r','LineWidth',3)
hold on
plot(centers,cumsum(dist_inhibited),'b','LineWidth',3)
plot(centers,cumsum(dist_tagged),'k','LineWidth',3)

%% pdf #2

tagged = find(inx_act&activation_start<3);
activated = setdiff(find(inx_act&~inx_inh),tagged);
inhibited = find(inx_inh&~inx_act);
ai = find(inx_act&inx_inh);
inx = activation_peak(ai) < inhibition_peak(ai);
activated_inhibited = ai(inx);
inhibited_activated = ai(~inx);

edges = 0:5:100;
edges = -2.5:5:102.5;
% edges = 0:3:100;
% edges = [0 1 2 4 8 16 32 64 128];
% edges = 0:10:100;
% edges = [0 3.^(0:5)];
centers = (edges(1:end-1) + edges(2:end)) / 2;
dist_activated = histc(activation_peak(activated),edges);
dist_activated = dist_activated(1:end-1);
dist_activated = dist_activated / sum(dist_activated);

dist_inhibited = histc(inhibition_peak(inhibited),edges);
dist_inhibited = dist_inhibited(1:end-1);
dist_inhibited = dist_inhibited / sum(dist_inhibited);

dist_activated_inhibited = histc(activation_peak(activated_inhibited),edges);
dist_activated_inhibited = dist_activated_inhibited(1:end-1);
dist_activated_inhibited = dist_activated_inhibited / sum(dist_activated_inhibited);

dist_inhibited_activated = histc(inhibition_peak(inhibited_activated),edges);
dist_inhibited_activated = dist_inhibited_activated(1:end-1);
dist_inhibited_activated = dist_inhibited_activated / sum(dist_inhibited_activated);

dist_tagged = histc(activation_peak(tagged),edges);
dist_tagged = dist_tagged(1:end-1);
dist_tagged = dist_tagged / sum(dist_tagged);

figure
plot(centers,dist_activated,'Color','r','LineWidth',3)
hold on
plot(centers,dist_inhibited,'Color','b','LineWidth',3)
plot(centers,dist_tagged,'Color','k','LineWidth',3)
plot(centers,dist_activated_inhibited,'Color',[0.6 0 0],'LineWidth',3)
plot(centers,dist_inhibited_activated,'Color',[0 0 0.6],'LineWidth',3)
plot(centers,dist_tagged,'Color','k','LineWidth',3)

%% cdf #2

figure
plot(centers,cumsum(dist_activated),'r','LineWidth',3)
hold on
plot(centers,cumsum(dist_inhibited),'b','LineWidth',3)
plot(centers,cumsum(dist_activated_inhibited),'Color',[0.6 0 0],'LineWidth',3)
plot(centers,cumsum(dist_inhibited_activated),'Color',[0 0 0.6],'LineWidth',3)
plot(centers,cumsum(dist_tagged),'k','LineWidth',3)

%% pdf #3 - accepted version

tagged = find(inx_act&activation_start<3);
activated = setdiff(find(inx_act&~inx_inh),tagged);
inhibited = find(inx_inh&~inx_act);
ai = find(inx_act&inx_inh);
inx = activation_peak(ai) < inhibition_peak(ai);
activated_inhibited = ai(inx);
inhibited_activated = ai(~inx);
activated = [activated activated_inhibited];
inhibited = [inhibited inhibited_activated];

edges = 0:5:100;
edges = -1:2:101;
% edges = 0:3:100;
% edges = [0 1 2 4 8 16 32 64 128];
% edges = 0:10:100;
% edges = [0 3.^(0:5)];
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
plot(centers,dist_tagged,'Color','k','LineWidth',3)

%% cdf #3

% figure
% plot(centers,cumsum(dist_activated),'r','LineWidth',3)
% hold on
% plot(centers,cumsum(dist_inhibited),'b','LineWidth',3)
% plot(centers,cumsum(dist_tagged),'k','LineWidth',3)

% figure
% plot(sort(activation_peak(activated),'ascend'),cumsum(ones(1,length(activated)))/length(activated),'r','LineWidth',3)
% hold on
% plot(sort(inhibition_peak(inhibited),'ascend'),cumsum(ones(1,length(inhibited)))/length(inhibited),'b','LineWidth',3)
% plot(sort(activation_peak(tagged),'ascend'),cumsum(ones(1,length(tagged)))/length(tagged),'k','LineWidth',3)

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
plot(centers,cumsum(dist_activated),'r','LineWidth',3)
hold on
plot(centers,cumsum(dist_inhibited),'b','LineWidth',3)
plot(centers,cumsum(dist_tagged),'k','LineWidth',3)