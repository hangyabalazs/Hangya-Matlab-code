%% Load

fileindex = 3;
[stimuli, spikes, param] = load_data;
if isequal(fileindex,1)
    pnn = 0;
else
    pnn = param(fileindex-1).nseconds;
end
[lfp,time] = get_trace('DataFile',param(fileindex).datafile,'Range',...
    [0 param(fileindex).nseconds],'Triggers',param(fileindex).position_offset_time+0.0001,'Channel',2);
sr = param(fileindex).samplerate;
tzero = param(fileindex).position_offset_time;   % starting point in sec

%% Plot raw data

figure
plot(time,lfp)
vdisc = (cell2mat({stimuli.trigger}) - tzero);
vdisc(vdisc>time(end)) = [];
hold on
line([vdisc; vdisc],[ones(1,length(vdisc))*2; ones(1,length(vdisc))*(-6)],'Color','g')


%% Stimulus parameters

dtf = cell2mat({stimuli(:).datafile});
tps = {stimuli(:).type};
trs = cell2mat({stimuli(:).trigger});

%% Filter stimuli

lightstim = strcmp(tps,'ledpulse');
% lightstim2 = strcmp(tps,'tone');
% lightstim2 = logical([1 lightstim2(1:end-1)]);
% lightstim = lightstim & lightstim2;
fl = find(lightstim);
dur = arrayfun(@(k)stimuli(k).param.duration,fl);
lightstim = fl(dur==40);

%% Downsample

dsr = 1000;
cst = sr / dsr;
lfp2 = lfp(1:cst:end);

%% Plot raw data

vdisc = (cell2mat({stimuli.trigger}) - tzero) * dsr;
vdisc(vdisc>length(lfp2))=[];
figure;plot(lfp2)
hold on;
line([vdisc; vdisc],[ones(1,length(vdisc))*2; ones(1,length(vdisc))*(-6)],'Color','g')

%% STA

vdisc = (cell2mat({stimuli(lightstim).trigger}) - tzero) * dsr;
vdisc = round(vdisc);

wn = 4 * dsr;
st = stacall(vdisc,lfp2,dsr,wn);

%% Filter STA

mns = min(st(:,1:floor(wn/2)),[],2);
sinx = find(mns<-2.5);
st2 = st;
st2(sinx,:) = [];

figure
imagesc(st2)
figure
plot(mean(st2))

eft = min(st(:,ceil(wn/2):ceil(wn/2)+0.1*dsr),[],2);
figure
plot(mns,eft,'.')

%% Compare noise with noise + light

noisestim = strcmp(tps,'whitenoise');
noiselightstim = cellfun(@(k)isequal(k,[{'ledpulse'} {'whitenoise'}]),tps);

vdisc = (cell2mat({stimuli(noisestim).trigger}) - tzero) * dsr;
vdisc = round(vdisc);
wn = 4 * dsr;
stn = stacall(vdisc,lfp2,dsr,wn);

vdisc = (cell2mat({stimuli(noiselightstim).trigger}) - tzero) * dsr;
vdisc = round(vdisc);
pulse_latency = arrayfun(@(k)stimuli(k).param{2}.pulse_latency,find(noiselightstim)) / 1000 * dsr;
wn = 4 * dsr;
stnl = stacall(vdisc+pulse_latency,lfp2,dsr,wn);

figure
plot(mean(stn))
hold on
plot(mean(stnl),'r')

%% Wavelet

dsr2 = 200;
cst2 = sr / dsr2;
lfp3 = lfp(1:cst2:end);

vdisc = (cell2mat({stimuli(lightstim).trigger}) - tzero) * dsr2;
vdisc = round(vdisc);

[pow,phase,f] = eegwavelet(lfp3,dsr2);
[ersp erspa] = ers(vdisc,pow,dsr2);

figure
imagesc(log(ersp))
b_rescaleaxis('Y',f)
setappdata(gca,'scaley',f)
b_zoomset_for_wavelet

figure
I = log(ersp) - log(repmat(mean(ersp,2),1,size(ersp,2)));
imagesc(I)
b_rescaleaxis('Y',f)
setappdata(gca,'scaley',f)
b_zoomset_for_wavelet

%% Filter wavelet

erspa(sinx,:,:) = [];   % filter
ersp2 = squeeze(mean(erspa,1));
figure
imagesc(log(ersp2))
b_rescaleaxis('Y',f)
setappdata(gca,'scaley',f)
b_zoomset_for_wavelet

%% Spectrogram

dsr2 = 400;
cst2 = sr / dsr2;
lfp3 = lfp(1:cst2:end);

vdisc = (cell2mat({stimuli(lightstim).trigger}) - tzero) * dsr2;

f = 1:(dsr2/4);
pow = log(abs(spectrogram(lfp3,round(110-1000/length(lfp3)),100,f,dsr2)).^2);
[ersp erspa] = ers(round(vdisc/(length(lfp3)/size(pow,2))),pow,dsr2/(length(lfp3)/size(pow,2)));

figure
imagesc(1:size(ersp,2),f,ersp)

%% Spectrogram2 - STA subtracted

sta = mean(st);
vdisc = (cell2mat({stimuli(lightstim).trigger}) - tzero) * dsr;
vdisc = round(vdisc);
lfp2_ = lfp2;
wn2 = round((length(sta)-1)/2);
for t = 1:length(vdisc)
    lfp2_(vdisc(t)-wn2:vdisc(t)+wn2) = lfp2_(vdisc(t)-wn2:vdisc(t)+wn2) - sta';
end

dsr2 = 500;
cst2 = dsr / dsr2;
lfp3 = lfp2_(1:cst2:end);
vdisc = (cell2mat({stimuli(lightstim).trigger}) - tzero) * dsr2;

f = 1:(dsr2/4);
pow = log(abs(spectrogram(lfp3,round(110-1000/length(lfp3)),100,f,dsr2)).^2);
[ersp erspa] = ers(round(vdisc/(length(lfp3)/size(pow,2))),pow,dsr2/(length(lfp3)/size(pow,2)));

figure
imagesc(1:size(ersp,2),f,ersp)