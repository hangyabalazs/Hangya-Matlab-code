% Extract spikes correlated with light-evoked spikes.

%% Load spikes

cellid = 'n045_121217x_4.6';  % closest to light-evoked spikes

% Load spikes from Ntt file
Nttfn = cellid2fnames(cellid,'Ntt');
[all_spikes wv] = LoadTT_NeuralynxNT(Nttfn);
TIMEFACTOR = getpref('cellbase','timefactor');    % scaling factor to convert spike times into seconds
all_spikes = all_spikes * TIMEFACTOR;
spk = loadcb(cellid,'Spikes');

% Restrict to the stim. period
SE = loadcb(cellid,'Stimevents');
ps1 = find(~isnan(SE.ProtocolStart));
ps2 = SE.ProtocolStart(ps1(end));
spk = spk(spk>ps2);
[jn1 jn2 spkinx] = intersect(spk,all_spikes);
if ~isequal(jn1,spk)  % check if all files were imported
    error('light_corr_cell:SpikeTimeMismatch','Mismatch between saved spike times and Ntt time stamps.')
end

%% Find largest channel

mean_wv = squeeze(nanmean(wv(spkinx,:,:),1));   % mean waveform
amx = max(max(mean_wv));     % absolute maximum of mean waveforms
[mx my] = find(mean_wv==amx,1,'first');     % mx: largest channel

%% bootstrap

spike_num = length(all_spikes);
sample_num = 1000;
bst = 1000;
wave_spont = squeeze(wv(:,mx,:));
wave_evoked = squeeze(mean(wv(spkinx,mx,:),1));
for k = 1:bst
    subset_idx = randsample(1:spike_num,sample_num);
    
    sr = 32552; % DigiLynx sampling rate
    rng = round(0.00075*sr);    % number of data points in 750 us (default censored period of DigiLynx)
    
    for j = 1:sample_num
        
        % correlation
        pr = corrcoef(wave_spont(subset_idx(j),1:rng),wave_evoked(1:rng));
        R_tmp(j) = pr(1,2);
        
        % dot product
        D_tmp(j) = dot(wave_spont(subset_idx(j),1:rng),wave_evoked(1:rng));
    end
    
    R(k) = mean(R_tmp);
    D(k) = mean(D_tmp);
end


%% thresholding R

thr = prctile(R,99.9);
for k = 1:spike_num
    pr2 = corrcoef(wave_spont(k,1:rng),wave_evoked(1:rng));
    R2(k) = pr2(1,2);
end
clusterinx = R2 > thr;
TS = all_spikes(clusterinx);

%% based on all corr. values

thr = prctile(R2,95);
clusterinx = R2 > thr;
TS = all_spikes(clusterinx);

%% based on difference

R2 = sum(abs(wave_spont-repmat(wave_evoked',spike_num,1)),2);
thr = prctile(R2,5);
clusterinx = R2 < thr;
TS = all_spikes(clusterinx);

%% based on diff, all tetrodes

wave_spont1 = squeeze(wv(:,1,:));
wave_evoked1 = squeeze(mean(wv(spkinx,1,:),1));
wave_spont2 = squeeze(wv(:,2,:));
wave_evoked2 = squeeze(mean(wv(spkinx,2,:),1));
wave_spont3 = squeeze(wv(:,3,:));
wave_evoked3 = squeeze(mean(wv(spkinx,3,:),1));
wave_spont4 = squeeze(wv(:,4,:));
wave_evoked4 = squeeze(mean(wv(spkinx,4,:),1));
R2 = sum(abs(wave_spont1-repmat(wave_evoked1',spike_num,1)),2) + ...
    sum(abs(wave_spont2-repmat(wave_evoked2',spike_num,1)),2) + ...
    sum(abs(wave_spont3-repmat(wave_evoked3',spike_num,1)),2) + ...
    sum(abs(wave_spont4-repmat(wave_evoked4',spike_num,1)),2);
thr = prctile(R2,5);
clusterinx = R2 < thr;
TS = all_spikes(clusterinx);


%% write

TS = TS * 1e4;
save('TT4_9.mat','TS');