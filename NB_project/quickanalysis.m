%%

animalID = 'n046';
animalID2 = 'nb046';
sessionID = '130108a'

fullpth = [getpref('cellbase','datapath') '\' animalID '\' sessionID '\'];

%%

nlxcsc2mat2(fullpth,'Channels','Events')

%%

TE = solo2trialevents4_auditory_gonogo_pulsepal([fullpth 'data_@auditory_gonogo_pulsepal_balazs_' animalID2 '_' sessionID '.mat']);

%%
TE = solo2trialevents4_auditory_gonogo([fullpth 'data_@auditory_gonogo_balazs_' animalID2 '_' sessionID '.mat']);

%%

MakeTrialEvents2_gonogo(fullpth)

%%

[ev, ep] = defineEventsEpochs_gonogo;

%%

% addnewsessions

addnewcells

%%

cellids = findcell('rat',animalID,'session',sessionID)

%%
tic;
problem_behav_cellid = [];
for iC=1:length(cellids),
    cellid=cellids(iC);
    disp(cellid)
    try
         prealignSpikes(cellid,'FUNdefineEventsEpochs',@defineEventsEpochs_gonogo,'filetype','event','ifsave',1,'ifappend',0) % filetype='event' for behavior protocol
    catch
        disp('Error in prealignSpikes.');
        problem_behav_cellid=[problem_behav_cellid cellid];
    end
end
toc;

%% is predictive?

for k = 1:length(cellids)
    figure
%     viewcell2b(cellids(k),'TriggerName','StimulusOn','SortEvent','Stimulu
%     sOff','eventtype','behav','ShowEvents',{{'StimulusOff'}},'Partitions','#ResponseType','window',[-5 5])
    viewcell2b(cellids(k),'TriggerName','StimulusOn','SortEvent','TrialStart','eventtype','behav','ShowEvents',{{'StimulusOff'}},'Partitions','#ResponseType','window',[-5 5])
%     viewcell2b(cellids(k),'TriggerName','StimulusOn','SortEvent','StimulusOff','eventtype','behav','ShowEvents',{{'StimulusOff'}},'Partitions','#ResponseType','window',[-0.6 0.6])
end

%% Hit & FA

for k = 1:length(cellids)
    figure
%     viewcell2b(cellids(k),'TriggerName','LeftPortIn','SortEvent','StimulusOn','eventtype','behav','ShowEvents',{{'StimulusOn'}},'Partitions','#ResponseType','window',[-5 5])
%     viewcell2b(cellids(k),'TriggerName','DeliverFeedback','SortEvent','PseudoStimulusOn','eventtype','behav','ShowEvents',{{'LeftPortIn'}},'Partitions','#ResponseType','window',[-5 5])
        
%     viewcell2b(cellids(k),'TriggerName','LeftPortIn','SortEvent','StimulusOn','eventtype','behav','ShowEvents',{{'StimulusOn'}},'Partitions','#ResponseType','window',[-0.6 0.6],'isadaptive',2,'dt',0.001)
%     viewcell2b(cellids(k),'TriggerName','DeliverFeedback','SortEvent','PseudoStimulusOn','eventtype','behav','ShowEvents',{{'LeftPortIn'}},'Partitions','#ResponseType','window',[-3 3])
    viewcell2b(cellids(k),'TriggerName','DeliverFeedback','SortEvent','PseudoStimulusOn','eventtype','behav','ShowEvents',{{'LeftPortIn'}},'Partitions','#ResponseType','window',[-0.3 0.3],'isadaptive',2,'dt',0.001)   % for figures; adaptive only works if the non-smoothed psth is plotted (manually pspth=psth;)
    viewcell2b(cellids(k),'TriggerName','LeftPortIn','SortEvent','StimulusOn','eventtype','behav','ShowEvents',{{'LeftPortIn'}},'Partitions','#ResponseType','window',[-0.3 0.3],'isadaptive',2,'dt',0.001)   % for figures; adaptive only works if the non-smoothed psth is plotted (manually pspth=psth;)
    viewcell2b(cellids(k),'TriggerName','StimulusOn','SortEvent','PseudoStimulusOn','eventtype','behav','ShowEvents',{{'LeftPortIn'}},'Partitions','#ResponseType','window',[-0.3 0.3],'isadaptive',2,'dt',0.001)   % for figures; adaptive only works if the non-smoothed psth is plotted (manually pspth=psth;)
    viewcell2b(cellids(k),'TriggerName','DeliverFeedback','SortEvent','PseudoStimulusOn','eventtype','behav','ShowEvents',{{'LeftPortIn'}},'Partitions','#ResponseType','window',[-0.1 0.2],...
        'dt',0.001,'sigma',0.001)
end

%% does it depend on stim intensity?

for k = 1:length(cellids)
    figure
    viewcell2b(cellids(k),'TriggerName','StimulusOn','SortEvent','StimulusOff','eventtype','behav','ShowEvents',{{'StimulusOff'}},'Partitions','#StimulusDuration','window',[-5 5])
end

%% lickraster

for k = 1:length(cellids)
    figure
    viewcell2b(cellids(k),'TriggerName','LickIn','SortEvent','StimulusOff','eventtype','behav','ShowEvents',{{'StimulusOff'}},'Partitions','#ResponseType','window',[-5 5],'Num2Plot',500)
end


%%
problem_behav_cellid=[];
problem_stim_cellid=[];
allcells=listtag('cell');
for iC=1:length(cellids),
    cellid=cellids(iC);
    pathname=cellid2fnames(cellid,'Sess');
%     try
%         MakeTrialEvents2(pathname);
%         script_behavior_analysis
%         %pause
%     catch
%         problem_behav_cellid=[problem_behav_cellid cellid];
%         disp(cellid);
%     end
    try
        MakeStimEvents2(pathname,'BurstStartNttl',4)
    catch
        problem_stim_cellid=[problem_stim_cellid cellid];
    end
    try
        SE=load([pathname filesep 'StimEvents']);
        if isnan(SE.PulseOn),
            MakeStimEvents2(pathname,'BurstStartNttl',2)
        end
    catch
    end
end

%% Make New Stim Events

tic;
problem_stim_cellid = [];
for iC=1:length(cellids),
    cellid=cellids(iC);
    disp(cellid)
    try
         prealignSpikes(cellid,'FUNdefineEventsEpochs',@defineEventsEpochs_laserstim,'filetype','stim','ifsave',1,'ifappend',0) % filetype='event' for behavior protocol
    catch
        disp('Error in prealignSpikes.');
        problem_stim_cellid=[problem_stim_cellid cellid];
    end
end
toc;

%% View Light Triggered Raster PSTH
% eventtype stim
TrigEvent='BurstOn';
% TrigEvent='OmitPulse';
SEvent='BurstOff';
FNum=2;
win=[-0.2 2.5];
parts='all';
parts='#BurstNPulse';
% parts='#ProtocolID2';
dt=0.001;
sigma=0.001;
PSTHstd='on';
ShEvent={{'PulseOn','PulseOff','BurstOff'}};
ShEvColors=hsv(length(ShEvent{1}));
ShEvColors=mat2cell(ShEvColors,ones(size(ShEvColors,1),1),3);
for iCell=1:length(cellids),
    cellid=cellids(iCell);
%     cellid=allcells(iCell);
%     try
        figure
        viewcell2b(cellid,'TriggerName',TrigEvent,'SortEvent',SEvent,'ShowEvents',ShEvent,'ShowEventsColors',{ShEvColors},...
        'FigureNum',FNum,'eventtype','stim','window',win,'dt',dt,'sigma',sigma,'PSTHstd',PSTHstd,'Partitions',parts,...
        'EventMarkerWidth',0,'PlotZeroLine','off')
%         pause(1)
%     catch
%     end
end

%%

auditory_gonogo_psychplot2(animalID,sessionID)
auditory_gonogo_psychplot2(animalID,sessionID,[],(150:300))

%%


lighton_trial = find(TE.LightStimulation2==1);
lightoff_trial = find(TE.LightStimulation2==0);
auditory_gonogo_psychplot2(animalID,sessionID,[],lighton_trial)
auditory_gonogo_psychplot2(animalID,sessionID,[],lightoff_trial)

%%


lighton_trial = find(TE.LightStimulation3==1);
lightoff_trial = find(TE.LightStimulation3==0);
auditory_gonogo_psychplot2(animalID,sessionID,[],lighton_trial)
auditory_gonogo_psychplot2(animalID,sessionID,[],lightoff_trial)

%%


lighton_trial = find(TE.LightStimulation3==1);
lightoff_trial = find(TE.LightStimulation3==0);
lm1 = 600;
lm2 = 800;
auditory_gonogo_psychplot2(animalID,sessionID,[],lighton_trial(lighton_trial>lm1&lighton_trial<lm2))
auditory_gonogo_psychplot2(animalID,sessionID,[],lightoff_trial(lightoff_trial>lm1&lightoff_trial<lm2))

%%

cellid = cellids{3};

tsegs_spont = findSegs2(cellid,0,1,'segfilter','spont');
[selts_spont,seltsind_spont,selisi_spont] = extractSegSpikes(cellid,tsegs_spont);
[Data_spont,wave_spont]=get_mclust_waveforms4(cellid,selts_spont,[],'chans','all');

wave_spont = Data_spont.PeakAlignedWaveform;
figure
subplot(221)
mn1 = nanmean(squeeze(wave_spont(:,1,:)));
sd1 = nanstd(squeeze(wave_spont(:,1,:)));
plot(mn1,'r')
hold on
plot(mn1+sd1,'k')
plot(mn1-sd1,'k')

subplot(222)
mn2 = nanmean(squeeze(wave_spont(:,2,:)));
sd2 = nanstd(squeeze(wave_spont(:,2,:)));
plot(mn2,'r')
hold on
plot(mn2+sd2,'k')
plot(mn2-sd2,'k')

subplot(223)
mn3 = nanmean(squeeze(wave_spont(:,3,:)));
sd3 = nanstd(squeeze(wave_spont(:,3,:)));
plot(mn3,'r')
hold on
plot(mn3+sd3,'k')
plot(mn3-sd3,'k')

subplot(224)
mn4 = nanmean(squeeze(wave_spont(:,4,:)));
sd4 = nanstd(squeeze(wave_spont(:,4,:)));
plot(mn4,'r')
hold on
plot(mn4+sd4,'k')
plot(mn4-sd4,'k')

%%

cellid = cellids{3};

tsegs_stim = findSegs2(cellid,0.001, 0.006,'segfilter','stim');   % 0.0013; 0.0018 for the cholinergic
[selts_stim,seltsind_stim,selisi_stim] = extractSegSpikes(cellid,tsegs_stim);
[Data_stim,wave_stim]=get_mclust_waveforms4(cellid,selts_stim,[],'chans','all');

wave_stim = Data_stim.PeakAlignedWaveform;
figure
subplot(221)
mn1 = nanmean(squeeze(wave_stim(:,1,:)));
sd1 = nanstd(squeeze(wave_stim(:,1,:)));
plot(mn1,'r')
hold on
plot(mn1+sd1,'k')
plot(mn1-sd1,'k')

subplot(222)
mn2 = nanmean(squeeze(wave_stim(:,2,:)));
sd2 = nanstd(squeeze(wave_stim(:,2,:)));
plot(mn2,'r')
hold on
plot(mn2+sd2,'k')
plot(mn2-sd2,'k')

subplot(223)
mn3 = nanmean(squeeze(wave_stim(:,3,:)));
sd3 = nanstd(squeeze(wave_stim(:,3,:)));
plot(mn3,'r')
hold on
plot(mn3+sd3,'k')
plot(mn3-sd3,'k')

subplot(224)
mn4 = nanmean(squeeze(wave_stim(:,4,:)));
sd4 = nanstd(squeeze(wave_stim(:,4,:)));
plot(mn4,'r')
hold on
plot(mn4+sd4,'k')
plot(mn4-sd4,'k')

%%

P=max(wave_stim,[],3);
V=min(wave_stim,[],3);
P1=squeeze(P(:,1));
P3=squeeze(P(:,3));
V1=squeeze(V(:,1));
V3=squeeze(V(:,3));
A1=P1-V1;
A3=P3-V3;

%%

stimwave = squeeze(wave_stim(:,2,:));
spontwave = squeeze(wave_spont(:,2,:));

[coeff_stim, scores_stim] = princomp(nan2zero(stimwave));
stimpc1 = scores_stim(:,1);

[coeff_spont, scores_spont] = princomp(nan2zero(spontwave));
spontpc1 = scores_spont(:,1);

LB = prctile(spontpc1,5);
UB = prctile(spontpc1,95);

badinx = find(stimpc1<LB|stimpc1>UB);
badinx2 = seltsind_stim(badinx);
global BADINX
BADINX = badinx2;


%%

prealignSpikes({cellid},'FUNdefineEventsEpochs',@defineEventsEpochs_laserstim,'filetype','stim','ifsave',1,'ifappend',0)

figure
viewcell2b(cellid,'TriggerName',TrigEvent,'SortEvent',SEvent,'ShowEvents',ShEvent,'ShowEventsColors',{ShEvColors},...
    'FigureNum',FNum,'eventtype','stim','window',win,'dt',dt,'sigma',sigma,'PSTHstd',PSTHstd,'Partitions',parts,...
    'EventMarkerWidth',0,'PlotZeroLine','off')

%%

wave_stim2 = wave_stim;
wave_stim2(badinx,:,:) = [];
figure
subplot(221)
mn1 = nanmean(squeeze(wave_stim2(:,1,:)));
sd1 = nanstd(squeeze(wave_stim2(:,1,:)));
plot(mn1,'r')
hold on
plot(mn1+sd1,'k')
plot(mn1-sd1,'k')

subplot(222)
mn2 = nanmean(squeeze(wave_stim2(:,2,:)));
sd2 = nanstd(squeeze(wave_stim2(:,2,:)));
plot(mn2,'r')
hold on
plot(mn2+sd2,'k')
plot(mn2-sd2,'k')

subplot(223)
mn3 = nanmean(squeeze(wave_stim2(:,3,:)));
sd3 = nanstd(squeeze(wave_stim2(:,3,:)));
plot(mn3,'r')
hold on
plot(mn3+sd3,'k')
plot(mn3-sd3,'k')

subplot(224)
mn4 = nanmean(squeeze(wave_stim2(:,4,:)));
sd4 = nanstd(squeeze(wave_stim2(:,4,:)));
plot(mn4,'r')
hold on
plot(mn4+sd4,'k')
plot(mn4-sd4,'k')

%% spike waveform in a given window relative to an event

% cellid = 'n029_120202a_3.5';

% Evoked spikes
lim1 = 0.015;
lim2 = 0.022;

TE = loadcb(cellid,'TrialEvents');
valid_trials = ~isnan(TE.FalseAlarm);

tsegs_evoked = rel2abstimes(cellid,[lim1 lim2],'trial','LeftPortIn','valid_trials',valid_trials);   % convert period to epochs relative to pulses
selts_evoked = extractSegSpikes(cellid,tsegs_evoked);   % find putative stimualated spikes
wave_evoked = extractSpikeWaveforms(cellid,selts_evoked,'chans','all');  % get waveforms for the extracted spikes

% Spontaneous spikes
tsegs_spont = rel2abstimes(cellid,[-1.5,0],'trial','LeftPortIn','valid_trials',valid_trials);   % extract 2s periods before bursts
selts_spont = extractSegSpikes(cellid,tsegs_spont);     % extract spontaneous spikes
wave_spont = extractSpikeWaveforms(cellid,selts_spont,'chans','all');    % get waveforms for the extracted spikes

%%

weds = wave_evoked;
wsds = wave_spont;

% Plot light-evoked waveforms
out.H_evoked = figure('Position',[624 126 1092 852]);
H = set_subplots(2,2,0.05,0.05,'XTick',[],'XLim',[1 size(wsds,3)]);
for sp = 1:4
    hold(H(sp),'on')
    plot(H(sp),transpose(squeeze(weds(:,sp,:))))
end
title('Light-evoked spike shape')

% Plot spontaneous waveforms
out.H_spont = figure('Position',[624 126 1092 852]);
H = set_subplots(2,2,0.05,0.05,'XTick',[],'XLim',[1 size(wsds,3)]);
for sp = 1:4
    hold(H(sp),'on')
    plot(H(sp),transpose(squeeze(wsds(:,sp,:))))
end
title('Spontaneous spike shape')

%%

% Average waveforms
mean_spont = squeeze(nanmean(wave_spont,1));
mean_evoked = squeeze(nanmean(wave_evoked,1));

% Compare waveforms
out.H_compare = figure('Position',[624 126 1092 852]);
H = set_subplots(2,2,0.05,0.05,'XTick',[],'XLim',[1 size(wsds,3)]);
for sp = 1:4
    hold(H(sp),'on')
    plot(H(sp),transpose(squeeze(wsds(:,sp,:))),'Color',[0.9 0.9 0.9])
    plot(H(sp),transpose(mean_spont(sp,:)),'Color','k','LineWidth',6)
    plot(H(sp),transpose(mean_evoked(sp,:)),'Color',[0 153 255]/255,'LineWidth',2)
end
title('Compare spont. and light-evoked spike shape')
