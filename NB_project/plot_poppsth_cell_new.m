%% load

global DATAPATH
load([DATAPATH 'NB\responsesorter_new4\allPSTH.mat'])

%% PSTH stats

% PSTH statistics
activation_peak = [allstats.activation_peak];   % peak time of activation
activation_start = [allstats.activation_start];   % activation onset time
activation_end = [allstats.activation_end];   % activation offset time
maxvalue = [allstats.maxvalue];   % maximal firing rate
Wpa = [allstats.Wpa];   % Mann-Whitney test for significant activation
inhibition_peak = [allstats.inhibition_peak];   % peak time of inhibition
inhibition_start = [allstats.inhibition_start];   % inhibition onset time
inhibition_end = [allstats.inhibition_end];   % inhibition offset time
minvalue = [allstats.minvalue];   % minimal firing rate
Wpi = [allstats.Wpi];   % Mann-Whitney test for significant inhibition
baseline = [allstats.baseline];   % baseline firing rate

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
% inhibited = [inhibited'; inhibited_activated'];

%% groups

% Identified NB cells
selstr = ['"ChAT+"==1&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})'];
ChAT = selectcell(selstr);   % cell IDs for ChAT cells
selstr = ['"pChAT+"==1&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})'];
pChAT = selectcell(selstr);  % putative
selstr = ['"ChAT+"==0&"pChAT+"==0&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})'];
NT = selectcell(selstr);   % cell IDs for ChAT cells

[nms inxa ChATinx] = intersect(ChAT,tags);   % indices for ChAT cells
[nms inxa pChATinx] = intersect(pChAT,tags);   % indices for pChAT cells
[nms inxa NTinx] = intersect(NT,tags);   % indices for unidentified cells
ChATactinx = intersect(ChATinx,activated);   % indices for activated ChAT cells
pChATactinx = intersect(pChATinx,activated);   % indices for activated pChAT cells
allChATactinx = [ChATactinx pChATactinx];   % indices for activated ChAT and pChAT cells
NTactinx = intersect(NTinx,activated);  % indices for activated non-ChAT cells
NTinhinx = intersect(NTinx,inhibited);  % indices for inhibited non-ChAT cells

%%  population PSTHs for groups

% Images proportional to number of cells in the groups
figure
imshow(allpsth(allChATactinx,:));colormap hot
set(gca,'CLim',[0 8])

figure
imshow(allpsth(NTinhinx,:));colormap hot
set(gca,'CLim',[0 8])

figure
imshow(allpsth(NTactinx,:));colormap hot
set(gca,'CLim',[0 8])

%% average PSTH

figure
hold on;
baseline_act = mean(allpsth(NTactinx,1:200));
baseline_inh = mean(allpsth(NTinhinx,1:200));
baseline_ChAT = mean(allpsth(allChATactinx,1:200));
mn_act = mean(baseline_act);
mn_inh = mean(baseline_inh);
mn_ChAT = mean(baseline_ChAT);
errorshade(time,mean(allpsth(NTactinx,:))-mn_act,std(allpsth(NTactinx,:))/sqrt(size(allpsth(NTactinx,:),1)),...
    'LineColor',[0.7 0.7 0.7],'ShadeColor',[0.7 0.7 0.7])
errorshade(time,mean(allpsth(allChATactinx,:))-mn_ChAT,std(allpsth(allChATactinx,:))/sqrt(size(allpsth(allChATactinx,:),1)),...
    'LineColor',[51 204 51]/255,'ShadeColor',[51 204 51]/255)
errorshade(time,mean(allpsth(NTinhinx,:))-mn_inh,std(allpsth(NTinhinx,:))/sqrt(size(allpsth(NTinhinx,:),1)),...
    'LineColor',[1 0.6 0.78],'ShadeColor',[1 0.6 0.78])

line([0 0],[-4 20],'Color','k','LineStyle',':')
line([-200 600],[0 0],'Color','k')
set(gca,'box','off',...
    'FontSize',16,'TickDir','out','XLim',[-200 600],'Ylim',[-2.6 8])
xlabel('Time')
ylabel('Normalized firing rate')