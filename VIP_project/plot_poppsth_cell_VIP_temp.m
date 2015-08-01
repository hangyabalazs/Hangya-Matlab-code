%% load

global DATAPATH
load([DATAPATH 'VIP\responsesorter_Hit_newwin\allPSTH.mat'],'allstats')
Wpa1 = [allstats.Wpa];
load([DATAPATH 'VIP\responsesorter_FA_newwin\allPSTH.mat'],'allstats')
Wpa2 = [allstats.Wpa];
Wpa = nan(size(Wpa1));
Wpa(Wpa1<0.01&Wpa2<0.01) = 0;


%% load 

global DATAPATH
load([DATAPATH 'VIP\responsesorter_FA_newwin\allPSTH.mat'])


%% groups

selstr = '"VIP+"==1&"validity"==1';
VIP = selectcell(selstr);   % cellIDs of VIP+ cells
selstr = '"VIP+"~=1&"validity"==1';
NT = selectcell(selstr);   % cellIDs of VIP- cells

activation_peak = [allstats.activation_peak];   % peak time of activation
activation_start = [allstats.activation_start];
activation_end = [allstats.activation_end];
maxvalue = [allstats.maxvalue];
% Wpa = [allstats.Wpa];   % Mann-Whitney test for significant activation
inhibition_peak = [allstats.inhibition_peak];   % peak time of inhibition
inhibition_start = [allstats.inhibition_start];
inhibition_end = [allstats.inhibition_end];
minvalue = [allstats.minvalue];
Wpi = [allstats.Wpi];   % Mann-Whitney test for significant inhibition
baseline = [allstats.baseline];

% Groups of activated and inhibited cells
inx_act = Wpa < 0.01;   % significant activation
inx_inh = Wpi < 0.01;   % significant inhibition
vldty = getvalue('validity')';
activated = find(inx_act&~inx_inh&vldty);    % indices for activated cells
inhibited = find(inx_inh&~inx_act&vldty);    % indices for inhibited cells
ai = find(inx_act&inx_inh&vldty);   % indices for activated and inhibited cells
inx = activation_peak(ai) < inhibition_peak(ai);
activated_inhibited = sort(ai(inx));
inhibited_activated = ai(~inx);
activated = [activated'; activated_inhibited'];   % categorize based on first significant effect
inhibited = [inhibited'; inhibited_activated'];

[nms inxa VIPinx] = intersect(VIP,tags);   % indices for VIP cells
VIPactinx = intersect(VIPinx,activated);   % indices for activated VIP cells
NTactinx = setdiff(activated,VIPactinx);  % indices for unidentified cells
numVIP = length(VIPactinx);   % number of VIP cells
numNT = length(NTactinx);   % number of unidentified cells
VIP_PSTH = allpsth(VIPactinx,:);    % PSTHs of VIP cells
[~, ~, NTinx] = intersect(NT,tags);   % indices for VIP- cells

% Light effects
lightinhinx = [2 18 25 30 59 63 73 75 83 86 88];
lightinhinx2 = [2 18 30 59 63 75 83 86 88];  % inh-act not included
lightactinx = [33 68 80 84 89];
lightactinx2 = [33 68 84 89];  % act-inh not included

%% colors

brown = [0.32 0.19 0.19];
purple = [0.48 0.06 0.89];
grey1 = [0.7 0.7 0.7];
grey2 = [0.4 0.4 0.4];
green = [0 0.8 0];
red = [0.8 0 0];

%%  population PSTHs for groups

figure;imagesc(allpsth(NTactinx,:));colormap hot
set(gca,'CLim',[0 5])

figure;imagesc(allpsth(inhibited,:));colormap hot
set(gca,'CLim',[0 5])

% inhibited = find(inx_inh&~inx_act);  % inhibited only, w/o inhibited-activated
% figure;imagesc(allpsth(inhibited,:));colormap hot
% set(gca,'CLim',[0 5])

figure;imagesc(allpsth(VIPactinx,:));colormap hot
set(gca,'CLim',[0 5])

% Images proportional to number of cells in the groups
% figure;imshow(allpsth(VIPactinx,:));colormap hot
% set(gca,'CLim',[0 5])
% figure;imshow(allpsth(inhibited,:));colormap hot
% set(gca,'CLim',[0 5])
% figure;imshow(allpsth(NTactinx,:));colormap hot
% set(gca,'CLim',[0 5])

%% sort

PSTHs = allpsth(NTactinx,:);
NumPsths = size(PSTHs,1);

[mx mxinx] = max(PSTHs,[],2);
[srt srtinx] = sort(mx,'descend');
% figure
% imagesc(time,1:NumPsths,PSTHs(srtinx,:))

figure;imagesc(PSTHs(srtinx,:));colormap hot
set(gca,'CLim',[0 5])


%% average PSTH - sign. cells

figure
hold on;
errorshade(time,mean(allpsth(NTactinx,:)),std(allpsth(NTactinx,:))/sqrt(size(allpsth(NTactinx,:),1)),...
    'LineColor',[0.7 0.7 0.7],'ShadeColor',[0.7 0.7 0.7])
errorshade(time,mean(allpsth(VIPactinx,:)),std(allpsth(VIPactinx,:))/sqrt(size(allpsth(VIPactinx,:),1)),...
    'LineColor',[51 204 51]/255,'ShadeColor',[51 204 51]/255)
errorshade(time,mean(allpsth(inhibited,:)),std(allpsth(inhibited,:))/sqrt(size(allpsth(inhibited,:),1)),...
    'LineColor',[1 0.6 0.78],'ShadeColor',[1 0.6 0.78])

line([0 0],[-4 20],'Color','k','LineStyle',':')
line([-200 600],[0 0],'Color','k')
set(gca,'box','off',...
    'FontSize',16,'TickDir','out','XLim',[-200 600],'Ylim',[-2 3])
xlabel('Time')
ylabel('Normalized firing rate')

%% average PSTH #2 - all cells

figure
hold on;
nm = sum(~all(isnan(allpsth(NTinx,:)),2));
errorshade(time,nanmean(allpsth(NTinx,:)),nanstd(allpsth(NTinx,:))/sqrt(nm),...
    'LineColor',[0.7 0.7 0.7],'ShadeColor',[0.7 0.7 0.7])
errorshade(time,mean(allpsth(VIPinx,:)),std(allpsth(VIPinx,:))/sqrt(size(allpsth(VIPinx,:),1)),...
    'LineColor',[51 204 51]/255,'ShadeColor',[51 204 51]/255)

line([0 0],[-4 20],'Color','k','LineStyle',':')
line([-200 600],[0 0],'Color','k')
set(gca,'box','off',...
    'FontSize',16,'TickDir','out','XLim',[-200 600],'Ylim',[-1.5 4])
xlabel('Time')
ylabel('Normalized firing rate')

%% average PSTH - grouped based on light-effect

figure
hold on;
errorshade(time,mean(allpsth(lightactinx2,:)),std(allpsth(lightactinx2,:))/sqrt(size(allpsth(lightactinx2,:),1)),...
    'LineColor',[0.7 0.7 0.7],'ShadeColor',[0.7 0.7 0.7])
errorshade(time,mean(allpsth(VIPactinx,:)),std(allpsth(VIPactinx,:))/sqrt(size(allpsth(VIPactinx,:),1)),...
    'LineColor',[51 204 51]/255,'ShadeColor',[51 204 51]/255)
errorshade(time,mean(allpsth(lightinhinx2,:)),std(allpsth(lightinhinx2,:))/sqrt(size(allpsth(lightinhinx2,:),1)),...
    'LineColor',[1 0.6 0.78],'ShadeColor',[1 0.6 0.78])

line([0 0],[-4 20],'Color','k','LineStyle',':')
line([-200 600],[0 0],'Color','k')
set(gca,'box','off',...
    'FontSize',16,'TickDir','out','XLim',[-200 600],'Ylim',[-1.5 4])
xlabel('Time')
ylabel('Normalized firing rate')

%% scatter plot

figure
semilogy(activation_peak(activated),maxvalue(activated)./baseline(activated),'o','MarkerEdgeColor',grey1,'MarkerFaceColor',grey1,'MarkerSize',8)
hold on
finh = minvalue(inhibited) ./ baseline(inhibited);
finh(finh==0) = 0.002;
semilogy(inhibition_peak(inhibited),finh,'o','MarkerEdgeColor',grey2,'MarkerFaceColor',grey2,'MarkerSize',8)
% semilogy(activation_peak(tagged),maxvalue(tagged)./baseline(tagged),'o','MarkerEdgeColor',[0 0.7 0],'MarkerFaceColor',[0 0.7 0],'MarkerSize',8)
% for k = 1:length(activated)
%     line([activation_start(activated(k)) activation_end(activated(k))],...
%         [maxvalue(activated(k))./baseline(activated(k)) maxvalue(activated(k))./baseline(activated(k))],'Color',grey1)
% end
% for k = 1:length(inhibited)
%     line([inhibition_start(inhibited(k)) inhibition_end(inhibited(k))],...
%         [minvalue(inhibited(k))./baseline(inhibited(k)) minvalue(inhibited(k))./baseline(inhibited(k))],'Color',grey2)
% end
semilogy(activation_peak(VIPactinx),maxvalue(VIPactinx)./baseline(VIPactinx),'o','MarkerEdgeColor',green,'MarkerFaceColor',green,'MarkerSize',8)
% for k = 1:length(VIPactinx)
%     line([activation_start(VIPactinx(k)) activation_end(VIPactinx(k))],...
%         [maxvalue(VIPactinx(k))./baseline(VIPactinx(k)) maxvalue(VIPactinx(k))./baseline(VIPactinx(k))],'Color',green)
% end

%% scatter plot - non-sign. included

figure
semilogy(activation_peak(activated),maxvalue(activated)./baseline(activated),'o','MarkerEdgeColor',grey1,'MarkerFaceColor',grey1,'MarkerSize',8)
hold on
finh = minvalue(inhibited) ./ baseline(inhibited);
finh(finh==0) = 0.002;
semilogy(inhibition_peak(inhibited),finh,'o','MarkerEdgeColor',grey2,'MarkerFaceColor',grey2,'MarkerSize',8)
% semilogy(activation_peak(tagged),maxvalue(tagged)./baseline(tagged),'o','MarkerEdgeColor',[0 0.7 0],'MarkerFaceColor',[0 0.7 0],'MarkerSize',8)
% for k = 1:length(activated)
%     line([activation_start(activated(k)) activation_end(activated(k))],...
%         [maxvalue(activated(k))./baseline(activated(k)) maxvalue(activated(k))./baseline(activated(k))],'Color',grey1)
% end
% for k = 1:length(inhibited)
%     line([inhibition_start(inhibited(k)) inhibition_end(inhibited(k))],...
%         [minvalue(inhibited(k))./baseline(inhibited(k)) minvalue(inhibited(k))./baseline(inhibited(k))],'Color',grey2)
% end
semilogy(activation_peak(VIPinx),maxvalue(VIPinx)./baseline(VIPinx),'o','MarkerEdgeColor',green,'MarkerFaceColor',green,'MarkerSize',8)
% for k = 1:length(VIPactinx)
%     line([activation_start(VIPactinx(k)) activation_end(VIPactinx(k))],...
%         [maxvalue(VIPactinx(k))./baseline(VIPactinx(k)) maxvalue(VIPactinx(k))./baseline(VIPactinx(k))],'Color',green)
% end
semilogy(activation_peak(NTinx),maxvalue(NTinx)./baseline(NTinx),'o','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',4)


%% format for ungrouping in PowerPoint

xl = xlim;
yl = ylim;
xd = get(gco,'XData');
yd = get(gco,'YData');
inx = xd > xl(1) & xd < xl(2) & yd > yl(1) & yd < yl(2);
% inx = xd > xl(1) & xd < xl(2) & yd > yl(1) & yd < 3.5;
set(gco,'XData',xd(inx),'YData',yd(inx))
