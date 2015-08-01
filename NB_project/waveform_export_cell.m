%% 1
inx=4;
figure('Position',[624   474   448   504])
plot(transpose(squeeze(wsds(:,inx,:))),'Color',[0.4 0.4 0.4])
ylim([-2 4.5]*10000)
hold on
plot(nanmean(squeeze(wsds(:,inx,:)),1),'Color',[0 153 255]/255,'LineWidth',3)
axis off

%% 2
inx=4;
figure('Position',[624   474   448   504])
plot(transpose(squeeze(wsds(:,inx,:))),'Color',[0.4 0.4 0.4])
ylim([-2 4.5]*10000)
hold on
plot(nanmean(squeeze(wsds(:,inx,:)),1),'Color',[1 1 .4],'LineWidth',3)
axis off

%% 5
inx=4;
figure('Position',[624   474   448   504])
plot(transpose(squeeze(wsds(:,inx,:))),'Color',[0.4 0.4 0.4])
ylim([-2 4.5]*10000)
hold on
plot(nanmean(squeeze(wsds(:,inx,:)),1),'Color',[0 1 1],'LineWidth',3)
axis off

%% 3
inx=4;
figure('Position',[624   474   448   504])
plot(transpose(squeeze(wsds(:,inx,:))),'Color',[0.4 0.4 0.4])
ylim([-2 4.5]*10000)
hold on
plot(nanmean(squeeze(wsds(:,inx,:)),1),'Color',[1 .4 .4],'LineWidth',3)
axis off

%% 7
inx=4;
figure('Position',[624   474   448   504])
plot(transpose(squeeze(wsds(:,inx,:))),'Color',[0.4 0.4 0.4])
ylim([-2 4.5]*10000)
hold on
plot(nanmean(squeeze(wsds(:,inx,:)),1),'Color',[.6 1 .4],'LineWidth',3)
axis off

%% 4
inx=4;
figure('Position',[624   474   448   504])
plot(transpose(squeeze(wsds(:,inx,:))),'Color',[0.4 0.4 0.4])
ylim([-2 4.5]*10000)
hold on
plot(nanmean(squeeze(wsds(:,inx,:)),1),'Color',[1 .6 1],'LineWidth',3)
axis off

%%

get(gco,'Color')
ans =
     0     1     1
get(gco,'Color')
ans =
    1.0000    0.4000    0.4000
get(gco,'Color')
ans =
    0.6000    1.0000    0.4000
get(gco,'Color')
ans =
    1.0000    0.6000    1.0000

    
%%
    
axis on
set(gca,'Color',[.8 .8 .8])

%%  raster

lnb=findobj(allchild(gca),'type','line','color','k');
set(lnb,'Color','w')
axis off
set(lnb,'Color','w','LineWidth',3)
lnr=findobj(allchild(gca),'facecolor','r');
set(lnr,'FaceColor',[169 158 103]/255)

%% PSTH

figure
alignfilter = 'FalseAlarm==1';
alignevent = 'DeliverFeedback'
[psth, spsth, ~, ~, spt, stats] = ultimate_psth(cellid,'trial',...
alignevent,wn,...
'dt',0.001,'display',false,'sigma',0.02,'parts','all','isadaptive',2,...
'event_filter','custom','filterinput',alignfilter,'maxtrialno',Inf,...
'baselinewin',[-0.5 0],'testwin',[0 0.5],'relative_threshold',0.1);
plot(psth)
set(gca,'TickDir','out')

%%

ChAT = selectcell(['"ChAT+"==1&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})']);  % identified ChAT+ cells
ChAT_Hit_PSTH = cell2mat(getvalue('Hit_PSTH',ChAT));
ChAT_FA_PSTH = cell2mat(getvalue('FA_PSTH',ChAT));

[nChAT_Hit_PSTH1 nChAT_FA_PSTH1] = deal(nan(size(ChAT_Hit_PSTH)));
for k = 1:length(ChAT)
    MN = mean(ChAT_FA_PSTH(k,:));
    SD = std(ChAT_FA_PSTH(k,:));
    nChAT_Hit_PSTH1(k,:) = (ChAT_Hit_PSTH(k,:) - MN) / SD;
    nChAT_FA_PSTH1(k,:) = (ChAT_FA_PSTH(k,:) - MN) / SD;
end

%%

ChAT_Hit_PSTH_stats = nancell2struct(getvalue('Hit_PSTH_stats',ChAT));
ChAT_FA_PSTH_stats = nancell2struct(getvalue('FA_PSTH_stats',ChAT));
ChAT_Hit_baseline = [ChAT_Hit_PSTH_stats.baseline];
ChAT_FA_baseline = [ChAT_FA_PSTH_stats.baseline];

[nChAT_Hit_PSTH nChAT_FA_PSTH] = deal(nan(size(ChAT_Hit_PSTH)));
for k = 1:length(ChAT)
    MN = mean(ChAT_FA_PSTH(k,:));
    SD = std(ChAT_FA_PSTH(k,:));
    nChAT_Hit_PSTH(k,:) = log(ChAT_Hit_PSTH(k,:) / ChAT_Hit_baseline(k));
    nChAT_FA_PSTH(k,:) = log(ChAT_FA_PSTH(k,:) / ChAT_FA_baseline(k));
end

%%

ChAT_Hit_Wpa = [ChAT_Hit_PSTH_stats.Wpa];
ChAT_FA_Wpa = [ChAT_FA_PSTH_stats.Wpa];
ChAT_Hit_Wpi = [ChAT_Hit_PSTH_stats.Wpi];
ChAT_FA_Wpi = [ChAT_FA_PSTH_stats.Wpi];

%%

alignfilter = 'Hit==1';
alignevent = 'LeftWaterValveOn';
[psth, spsth, ~, ~, spt, stats] = ultimate_psth(ChAT(13),'trial',...
    alignevent,wn,...
    'dt',0.001,'display',true,'sigma',0.02,'parts','all','isadaptive',2,...
    'event_filter','custom','filterinput',alignfilter,'maxtrialno',Inf,...
    'baselinewin',[-0.5 0],'testwin',[0 0.5],'relative_threshold',0.1);

%%

figure
imagesc(ChAT_Hit_PSTH(2:5,:))
colormap(hot)
set(gca,'clim',[5 150])


figure
imagesc(ChAT_FA_PSTH([2:4 6:14],:))
colormap(hot)
set(gca,'clim',[5 150])


%%

figure
plot(mean(nChAT_Hit_PSTH1))
ylim([-1 10])

figure
plot(mean(nChAT_FA_PSTH1))
ylim([-1 10])

%%

figure
errorshade(-600:600,mean(nChAT_Hit_PSTH1),std(nChAT_Hit_PSTH1)/sqrt(size(nChAT_Hit_PSTH1,1)),...
    'LineColor',[51 204 51]/255,'ShadeColor',[51 204 51]/255)
ylim([-1 11])

figure
errorshade(-600:600,mean(nChAT_FA_PSTH1),std(nChAT_FA_PSTH1)/sqrt(size(nChAT_FA_PSTH1,1)),...
    'LineColor',[216 41 0]/255,'ShadeColor',[216 41 0]/255)
ylim([-1 11])

%% attention cell

lnb = findobj(allchild(gca),'type','line','color','k');
set(lnb,'Color','w')
lnr=findobj(allchild(gca),'facecolor','r');
set(lnr,'FaceColor',[169 158 103]/255,'EdgeColor',[169 158 103]/255)
newch = [lnr; allchild(gca)];
fn = find(newch==lnr);
newch(fn(2)) = [];
set(gca,'Children',newch)

%% more attention cells

set(gcf,'Position',[198 457 1101 498])
line([0 0],ylim,'Color',[0.9255    0.8392    0.8392],'LineWidth',3)
set(gca,'XColor','w','YColor','w','FontSize',16,'TickDir','out','Box','off','XTick',[-2 -1 0 1 2])
xlabel('Time - StimulusOn','FontSize',16,'Color','w')
ylabel('Rate (Hz)','FontSize',16,'Color','w')