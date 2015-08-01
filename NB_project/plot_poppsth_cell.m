%% load

global DATAPATH
load([DATAPATH 'NB\responsesorter3\allPSTH.mat'])

%% groups

selstr = ['"ChAT+"==1&"ID_PC">20&"Lr_PC"<0.15&' ...
    'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})'];
ChAT = selectcell(selstr);   % cellIDs of ChAT+ cells
pChAT = {'n029_120210a_3.3' ...
    'n029_120215a_3.4' 'n029_120221b_6.1' 'n029_120222b_4.1'};   % cellIDs of putative ChAT+ cells
I = [ChAT pChAT];
allChAT = I;
selstr = ['"PV+"==1&"ID_PC">20&"Lr_PC"<0.15&' ...
    'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})'];
PV = selectcell(selstr);   % cellIDs of PV+ cells

activation_peak = [allstats.activation_peak];   % peak time of activation
Wpa = [allstats.Wpa];   % Mann-Whitney test for significant activation
inhibition_peak = [allstats.inhibition_peak];   % peak time of inhibition
Wpi = [allstats.Wpi];   % Mann-Whitney test for significant inhibition

% Groups of activated and inhibited cells
inx_act = Wpa < 0.01;   % significant activation
inx_inh = Wpi < 0.01;   % significant inhibition
activated = find(inx_act&~inx_inh);    % indices for activated cells
inhibited = find(inx_inh&~inx_act);    % indices for inhibited cells
ai = find(inx_act&inx_inh);
inx = activation_peak(ai) < inhibition_peak(ai);
activated_inhibited = sort(ai(inx));
inhibited_activated = ai(~inx);
activated = [activated'; activated_inhibited'];   % categorize based on first significant effect
inhibited = [inhibited'; inhibited_activated'];

[nms inxa PVinx] = intersect(PV,tags);   % indices for PV cells
[nms inxa ChATinx] = intersect(ChAT,tags);   % indices for ChAT cells
PVactinx = intersect(PVinx,activated);   % indices for activated PV cells
ChATactinx = intersect(ChATinx,activated);   % indices for activated ChAT cells
NTactinx = setdiff(activated,union(ChATactinx,PVactinx));  % indices for unidentified cells
numChAT = length(ChATactinx);   % number of ChAT cells
numPV = length(PVactinx);   % number of PV cells
numNT = length(NTactinx);   % number of unidentified cells
ChAT_PSTH = allpsth(ChATactinx,:);   % PSTHs of ChAT cells
PV_PSTH = allpsth(PVactinx,:);    % PSTHs of PV cells
pChATactinx = [259 263 427 384 428 461 464 543 559];

activatedinx = setdiff(activated,pChATactinx);

%%  population PSTHs fo groups

figure;imagesc(allpsth(activatedinx,:));colormap hot
set(gca,'CLim',[0 5])

figure;imagesc(allpsth(inhibited,:));colormap hot
set(gca,'CLim',[0 5])

inhibited = find(inx_inh&~inx_act);  % inhibited only, w/o inhibited-activated
figure;imagesc(allpsth(inhibited,:));colormap hot
set(gca,'CLim',[0 5])

figure;imagesc(allpsth(pChATactinx,:));colormap hot
set(gca,'CLim',[0 5])

% Images proportional to number of cells in the groups
figure;imshow(allpsth(allChATactinx,:));colormap hot
set(gca,'CLim',[0 8])
figure;imshow(allpsth(inhibited,:));colormap hot
set(gca,'CLim',[0 8])
figure;imshow(allpsth(activatedinx,:));colormap hot
set(gca,'CLim',[0 8])

%% labels for enlarged image of allChAT population PSTH

set(gca,'FontSize',16)
xlabel('Time from response onset (s)','Color','w','FontSize',16)
ylabel('Neuron#','Color','w','FontSize',16)

%% average PSTH

figure
hold on;
% errorshade(time,mean(allpsth(PVactinx,:)),std(allpsth(PVactinx,:))/sqrt(size(allpsth(PVactinx,:),1)),...
%     'LineColor',[216 41 0]/255,'ShadeColor',[216 41 0]/255)
baseline_act = mean(allpsth(activatedinx,1:200));
baseline_inh = mean(allpsth(inhibited,1:200));
baseline_ChAT = mean(allpsth(allChATactinx,1:200));
mn_act = mean(baseline_act);
mn_inh = mean(baseline_inh);
mn_ChAT = mean(baseline_ChAT);
errorshade(time,mean(allpsth(activatedinx,:))-mn_act,std(allpsth(activatedinx,:))/sqrt(size(allpsth(activatedinx,:),1)),...
    'LineColor',[0.7 0.7 0.7],'ShadeColor',[0.7 0.7 0.7])
errorshade(time,mean(allpsth(allChATactinx,:))-mn_ChAT,std(allpsth(allChATactinx,:))/sqrt(size(allpsth(allChATactinx,:),1)),...
    'LineColor',[51 204 51]/255,'ShadeColor',[51 204 51]/255)
errorshade(time,mean(allpsth(inhibited,:))-mn_inh,std(allpsth(inhibited,:))/sqrt(size(allpsth(inhibited,:),1)),...
    'LineColor',[1 0.6 0.78],'ShadeColor',[1 0.6 0.78])

line([0 0],[-4 20],'Color','k','LineStyle',':')
line([-200 600],[0 0],'Color','k')
set(gca,'box','off',...
    'FontSize',16,'TickDir','out','XLim',[-200 600],'Ylim',[-2.6 8])
xlabel('Time')
ylabel('Normalized firing rate')

%% format for ungrouping in PowerPoint

xl = xlim;
yl = ylim;
xd = get(gco,'XData');
yd = get(gco,'YData');
inx = xd > xl(1) & xd < xl(2) & yd > yl(1) & yd < yl(2);
% inx = xd > xl(1) & xd < xl(2) & yd > yl(1) & yd < 3.5;
set(gco,'XData',xd(inx),'YData',yd(inx))
