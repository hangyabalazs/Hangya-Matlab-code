%% spike waveform in a given window relative to an event

cellid = 'n029_120202a_3.5';

% Evoked spikes
lim1 = 0.015;
lim2 = 0.018;

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
