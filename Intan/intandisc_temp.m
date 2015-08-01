function intandisc(data,thr)

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
    cvdisc = cell(1,4);
    for iC = 1:4
        chnum = (iT - 1) * 4 + iC;  % current channel
        disp(['tetrode: ' num2str(iT) '   Channel: ' num2str(chnum)])
        cdata = data(:,chnum)';   % data from one channel
        cvdisc{iC} = disc(-cdata,thr);   % discriminate (spike times)
    end
    tvdisc{iT} = cell2mat(cvdisc);
    tvdisc{iT} = sort(tvdisc{iT},'ascend');  % all spike times from one tetrode
        
    % Dead time
    dtv = diff(tvdisc{iT});   % ISI
    censor_inx = dtv < dtp;   % ISI < dead time
    tvdisc{iT}(find(censor_inx)+1) = [];   % remove spikes within the censor period
    
    % Waveform
    winx = repmat(tvdisc{iT}(:)+winp(1),1,sum(abs(winp))) + repmat(0:sum(abs(winp))-1,length(tvdisc{iT}),1);
    dtmp = data(:,(iT-1)*4+1);
    wv1 = dtmp(winx);
    dtmp = data(:,(iT-1)*4+2);
    wv2 = dtmp(winx);
    dtmp = data(:,(iT-1)*4+3);
    wv3 = dtmp(winx);
    dtmp = data(:,(iT-1)*4+4);
    wv4 = dtmp(winx);
    1;
    
    
end

keyboard