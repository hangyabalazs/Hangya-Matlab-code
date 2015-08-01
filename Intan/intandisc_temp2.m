function intandisc(data,thr)

% Results directory
global DATADIR
resdir = [DATADIR 'Intan' filesep 'nba47' filesep];
savestr0 = 'save TT';

% Sampling rate
sr = 25000;   % 25 kHz
deadtime = 0.0005;   % 500 us dead time
dtp = deadtime * sr;   % dead time in data points

% Waveform window
win = [-0.00024 0.00096];   % -240 to 960 us
winp = win * sr;   % waveform window in data points

% Threshold discrimintaion
NumChannels = size(data,2);   % number of channels - 16
NumTetrodes = NumChannels / 4;   % number of tetrodes - 4
tvdisc = cell(1,NumTetrodes);
for iT = 1:NumTetrodes
    [cvdisc tdata] = deal(cell(1,4));
    for iC = 1:4
        chnum = (iT - 1) * 4 + iC;  % current channel
        disp(['tetrode: ' num2str(iT) '   Channel: ' num2str(chnum)])
        cdata = data(:,chnum)';   % data from one channel
        cvdisc{iC} = disc(-cdata,thr);   % discriminate (spike times)
        tdata{iC} = cdata;   % data from the current tetrode
    end
    tvdisc{iT} = cell2mat(cvdisc);
    tvdisc{iT} = sort(tvdisc{iT},'ascend');  % all spike times from one tetrode
        
    % Dead time
    dtv = diff(tvdisc{iT});   % ISI
    censor_inx = dtv < dtp;   % ISI < dead time
    tvdisc{iT}(find(censor_inx)+1) = [];   % remove spikes within the censor period
    
    % Waveform
    winx = repmat(tvdisc{iT}(:)+winp(1),1,sum(abs(winp))) + repmat(0:sum(abs(winp))-1,length(tvdisc{iT}),1);
    wv = nan(size(winx,1),4,size(winx,2));   % waveforms: cells x channels x time
    for iC = 1:4
        wv(:,iC,:) = tdata{iC}(winx);   % waveform data
    end
    eval(['WaveForm' num2str(iT) ' = wv;']);
%     wv1 = tdata{1}(winx);
%     wv2 = tdata{2}(winx);
%     wv3 = tdata{3}(winx);
%     wv4 = tdata{4}(winx);
    
    % Spike times
    spike_times = tvdisc{iT} / sr; %#ok<NASGU>
    eval(['TimeStamps' num2str(iT) ' = spike_times;']);
    savestr = [savestr0 num2str(iT) ' WaveForm' num2str(iT) ' TimeStamps' num2str(iT)];
    
    % Save
    cdr = pwd;
    cd(resdir)
    eval(savestr)
    cd(cdr)
end

keyboard