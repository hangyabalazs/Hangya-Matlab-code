function snwphase4
%SNWPHASE4   Whisking phase relative to respiration.
%   SNWPHASE4 calculates phase distributions of whisking events (EMG)
%   relative to respiratory cycles.
%
%   Raw data is resampled at 1 kHz. Respiration (sniffing) is low-pass
%   filtered at 30 Hz. Respiration phase is calculated by
%   Hilbert-transform. Whisking is high-pass filtered at 10 Hz, smoothed by
%   10 ms window root-mean-square, high-pass filtered at 50 Hz and
%   rectified. Individual whisks are detected by threshold disrimintaion.
%   Different whisking patterns are separated (number of whisks:sniff
%   cycles, 3:1, 2:1, 1:1, 1:2, 1:3) and phase histograms restricted to the
%   patterns are plotted and saved.
%
%   See also SSPW3 and PHASEDEP.

% Define results directory
global DATAPATH
resdir = [DATAPATH 'SniffWhisk\phasehist2_\'];

% Data import
plotting = true;
global DATADIR
dsr = 1000;   % downsample rate
rrat = 'R7';
rdate = '2008-09-03_16-51-55';
csc1 = 'CSC4.ncs';
csc2 = 'CSC7.ncs';
msgn = -1;
wl = dsr / 100;     % 10 ms windows for RMS (for 1 kHz samp. rate)
% rrat = 'P5';
% rdate = 'P5 2005-11-8';
% csc1 = 'CSC4.ncs';
% csc2 = 'CSC7.ncs';
% msgn = 1;
% wl = dsr / 50;     % 10 ms windows RMS (for 1 kHz samp. rate) (/50 for 20 ms)

fn = ([DATADIR 'HSW\' rrat '\' rdate '\']);
[TimeStamp, ChanNum, SampleFrequency, NumValSamples, Samples, NlxHeader] = ...
    Nlx2MatCSC([fn csc1],[1 1 1 1 1],1,1,1);  % sniffing
sniff = msgn * Samples(:);
sr_sniff = 512 / mean(diff(TimeStamp)) * 10^6;    % original sampling rate (changed later)
dt_sniff = 1 / sr_sniff;
starttime1 = TimeStamp(1) / 10^6;
time_sniff = (0:length(sniff)-1) * dt_sniff + starttime1;
[TimeStamp, ChanNum, SampleFrequency, NumValSamples, Samples, NlxHeader] = ...
    Nlx2MatCSC([fn csc2],[1 1 1 1 1],1,1,1);  % whisking
whisk = Samples(:);
sr_whisk = 512 / mean(diff(TimeStamp)) * 10^6;    % original sampling rate for EMG (different from that for respiration in R7)
dt_whisk = 1 / sr_whisk;
starttime2 = TimeStamp(1) / 10^6;
time_whisk = (0:length(whisk)-1) * dt_whisk + starttime2;
timediff = starttime2 - starttime1;
% sniff = linterp(time_sniff,sniff,time_whisk);

% Resample
[p q] = rat(dsr/sr_sniff);    % resample respiration at 1000 Hz
sniff = resample(sniff,p,q);
[p q] = rat(dsr/sr_whisk);    % resample EMG at 1000 Hz
whisk = resample(whisk,p,q);
sr = dsr;  % new sampling rate
if timediff > 0.0001       % deal with misaligned data
    if ~isequal(sr,1000)
        error('Data realignment only works for samplin rate = 1000 Hz.')
    end
    leadtime = timediff * sr;    % lead of sniff in ms
    sniff(1:round(leadtime)) = [];  % drop lead: now they co-start
elseif timediff < -0.0001
    if ~isequal(sr,1000)
        error('Data realignment only works for samplin rate = 1000 Hz.')
    end
    lagtime = (-1) * timediff * sr;    % lag of sniff in ms
    whisk(1:round(lagtime)) = [];  % drop lead of whisk: now they co-start
end
ln = min(length(sniff),length(whisk));
sniff = sniff(1:ln);
whisk = whisk(1:ln);   % now they co-terminate

% Deal with individual artifacts
if isequal(rdate,'2008-09-03_16-42-14')
    whisk = whisk(1:2.25*10^5);
    sniff = sniff(1:2.25*10^5);
elseif isequal(rdate,'2008-09-21_21-08-00')
    whisk = whisk(1:1.54*10^6);
    sniff = sniff(1:1.54*10^6);
end

% Filter
nqf = sr / 2;
flt = fir1(1024,30/nqf,'low');
fsniff = filtfilt(flt,1,sniff);     % lowpass filter sniffing at 30 Hz
flt = fir1(4096,10/nqf,'high');
fwhisk = filtfilt(flt,1,whisk);     % highpass filter whisking at 10 Hz

% Respiration phase
resp_phase = angle(hilbert(standardize(fsniff)));    % Hilbert-transform

% Downsampled resp. phase
lenw = length(fwhisk);
lle = floor(lenw/wl) * wl;
resp_phase2 = resp_phase(1:lle);
resp_phase2 = resp_phase2(round(wl/2):wl:end);   % downsampled resp. phase aligned to whisking RMS

% Detect whisking
% [st nd pk fwhisk_ rms_whisk thr] = sspw3(whisk,sr,10,10,inf,2,1);
% [st nd pk fwhisk_ rms_whisk thr] = sspw3(whisk,sr,wl,10,inf,[],[]);  % whisking detection, see SSPW3
[st nd pk0 pk fwhisk_ rms_whisk thr] = sspw3_temp(whisk,sr,wl,10,inf,[],[],500);  % whisking detection, see SSPW3
if plotting
    figure
    plot(rms_whisk)
    hold on
    line([pk/wl; pk/wl],[zeros(size(pk)); ones(size(pk))*10000],'Color','green')
end

% Get resp. cycles
fn0 = valuecrossing(1:length(resp_phase2),resp_phase2',0,'down');   % limits for resp. cycles
fn = round(fn0);
if plotting
    figure
    line([fn; fn],[ones(size(fn))-4; ones(size(fn))+3],'Color','green')
    hold on
    plot(resp_phase2)
end

% Whisking events
vdisc = round(pk/wl);
if plotting
    figure
    line([vdisc; vdisc],[zeros(size(pk)); ones(size(pk))*2000],'Color','green')
    hold on
    plot(fsniff(1:wl:end))
end

% Sort whisking events
grp = zeros(size(vdisc));     % codes for types of whisking events 
time_grp = zeros(size(resp_phase2));   % codes time-wise 
time_grp2 = zeros(size(resp_phase2));   % indicator of cycles (time) before 1:2 whisking events
for k = 8:length(vdisc)-12
    ffv = find(fn<vdisc(k),1,'last');    % fn: limits of respiration cycles
    fn_4pre = fn(ffv-3);     % cycle limit 4 before current events 
    fn_pppre = fn(ffv-2);    % cycle limit 3 before current events
    fn_ppre = fn(ffv-1);     % cycle limit 2 before current events
    fn_pre = fn(ffv);        % cycle limit 1 before current events
    fn_post = fn(ffv+1);     % cycle limit 1 after current events
    fn_ppost = fn(ffv+2);    % cycle limit 2 after current events
    fn_pppost = fn(ffv+3);   % cycle limit 3 after current events
    fn_4post = fn(ffv+4);    % cycle limit 4 after current events
    lv_pppre = length(vdisc(vdisc>fn_4pre&vdisc<fn_pppre));    % whisking events in 3 cycles before
    lv_ppre = length(vdisc(vdisc>fn_pppre&vdisc<fn_ppre));     % whisking events in 2 cycles before
    lv_pre = length(vdisc(vdisc>fn_ppre&vdisc<fn_pre));        % whisking events in 1 cycle before
    lv = length(vdisc(vdisc>fn_pre&vdisc<fn_post));            % whisking events in current cycle
    lv_post = length(vdisc(vdisc>fn_post&vdisc<fn_ppost));     % whisking events in 1 cycle after
    lv_ppost = length(vdisc(vdisc>fn_ppost&vdisc<fn_pppost));  % whisking events in 2 cycles after
    lv_pppost = length(vdisc(vdisc>fn_pppost&vdisc<fn_4post)); % whisking events in 3 cycles after
    switch lv
        case 2
            grp(k) = 2;     % 2:1 coupling
            time_grp(floor(fn_pre):ceil(fn_post)) = 2;
        case 3              % 3:1 coupling
            grp(k) = 3;
            time_grp(floor(fn_pre):ceil(fn_post)) = 3;
        case 1              % 1:1 coupling
            if lv_pre >= 1 || lv_post >= 1   % at least one event in 1 cycle pre or post
                grp(k) = (1);
                time_grp(floor(fn_pre):ceil(fn_post)) = 1;
            else            % 1:2 coupling
                if lv_ppre >= 1 || lv_ppost >= 1    % no event in pre and post, but at least one event 2 cycles away
                    grp(k) = 0.5;
                    time_grp(floor(fn_ppre):ceil(fn_ppost)) = 0.5;
                    time_grp2(floor(fn_ppre):ceil(fn_pre-1)) = 1;   % indicator for cyles before and after the ones with the whisking events in 1:2 coupling
                    time_grp2(floor(fn_post):ceil(fn_ppost-1)) = 1;
                else        % 1:3 coupling
                    if lv_pppre >= 1 || lv_pppost >= 1   % no event in 2 cycles around, but at least one event 3 cycles away
                        grp(k) = 1/3;
                        time_grp(floor(fn_pppre):ceil(fn_pppost)) = 1/3;
                    end
                end
            end
    end
end

% Eliminate short segments
dfs = find(diff(time_grp));
for k = 1:length(dfs)-1     % loop through uniform segments
    inx1 = dfs(k) + 1;  % segment start
    inx2 = dfs(k+1);    % segment end
    cseg = time_grp(inx1:inx2);   % segment (time domain)
    tcseg = length(cseg) / sr * wl;   % length of the segment in seconds
    cgrp = grp(vdisc>=inx1&vdisc<=inx2&grp~=0);   % whisking events within the segment
    if tcseg < 0.5 || length(cgrp) < 3  % exculde segments < 0.5 s and segments with < 3 whisking events
        time_grp(inx1:inx2) = 0;
        time_grp2(inx1:inx2) = 0;
        grp(vdisc>=inx1&vdisc<=inx2) = 0;
    end
end

% Split 1:1 coupling
grp_split = grp;   % mark every odd whisk in 1:1 coupling bouts
time_grp_split = time_grp;   % same in time
for k = 1:length(grp)-1
    ffv = find(fn<vdisc(k),1,'last');
    fn_pre = fn(ffv);        % cycle limit 1 before current events
    fn_post = fn(ffv+1);     % cycle limit 1 after current events
    if grp(k) == 1 && (grp(k-1) ~= 1 || grp_split(k-1) == 1)
        grp_split(k) = 0.9;   % code for odd whisk events in 1:1 coupling bouts
        time_grp_split(floor(fn_pre):ceil(fn_post)) = 0.9;
    end
end

% Respiration phase - whisking relationship
keyboard
time_grp_orig = time_grp;
rms_whisk_orig = rms_whisk;
rms_whisk(rms_whisk<0) = 0;     % rectify RMS: due to filtering, below-0 values could appear
time_grp1 = time_grp;   % 1:1 coupling
time_grp1(time_grp~=1) = 0;     % indicator variable for 1:1 coupling (time)
time_grp1(time_grp==1) = 1;
[R p] = dphist(time_grp1,rms_whisk,resp_phase2);    % whisking RMS conditioned on resp. phase
fnm = [resdir 'PHASE1TO1_' rrat '_' rdate '.fig'];   % save figure
saveas(gcf,fnm)
fnm = [resdir 'PHASE1TO1_' rrat '_' rdate '.mat'];   % save variables
save(fnm,'R','p','time_grp1')

time_grp1 = time_grp;   % 1:2 coupling
time_grp1(time_grp~=0.5) = 0;     % indicator variable for 1:2 coupling (time)
time_grp1(time_grp==0.5) = 1;
time_grp1(time_grp2==1) = 0;    % disregard those not having whisking
[R p] = dphist2(time_grp1,time_grp2,rms_whisk,resp_phase2);    % whisking RMS conditioned on resp. phase
fnm = [resdir 'PHASE1TO2_' rrat '_' rdate '.fig'];   % save figure
saveas(gcf,fnm)     % plot shows every cycle before/after the detected events in black, see 'time_grp2'
fnm = [resdir 'PHASE1TO2_' rrat '_' rdate '.mat'];   % save variables
save(fnm,'R','p','time_grp1','time_grp2')

time_grp1 = time_grp;   % 2:1 coupling
time_grp1(time_grp~=2) = 0;     % indicator variable for 2:1 coupling (time)
time_grp1(time_grp==2) = 1;
[R p] = dphist(time_grp1,rms_whisk,resp_phase2);    % whisking RMS conditioned on resp. phase
fnm = [resdir 'PHASE2TO1_' rrat '_' rdate '.fig'];   % save figure
saveas(gcf,fnm)     % plot shows every cycle before the detected events in black, see 'time_grp2'
fnm = [resdir 'PHASE2TO1_' rrat '_' rdate '.mat'];   % save variables
save(fnm,'R','p','time_grp1')

% 1:1 couling, split into two groups to compare with 1:2
time_grp_split1 = time_grp_split;   % 1:1 coupling, split
time_grp_split2 = time_grp_split;
time_grp_split1(time_grp_split~=0.9) = 0;     % indicator variable for odd cycles (time)
time_grp_split1(time_grp_split==0.9) = 1;
time_grp_split2(time_grp_split~=1) = 0;     % indicator variable for even cycles (time)
time_grp_split2(time_grp_split==1) = 1;
[R p] = dphist2(time_grp_split1,time_grp_split2,rms_whisk,resp_phase2);    % whisking RMS conditioned on resp. phase
fnm = [resdir 'PHASE1TO1SPLIT_' rrat '_' rdate '.fig'];   % save figure
saveas(gcf,fnm)     % plot shows every even cycle in black
fnm = [resdir 'PHASE1TO1SPLIT_' rrat '_' rdate '.mat'];   % save variables
save(fnm,'R','p','time_grp_split1','time_grp_split2')

% Save additional variables
fnm = [resdir 'PHASE_' rrat '_' rdate '.mat'];   % save variables
save(fnm,'rms_whisk','rms_whisk_orig','vdisc','thr','resp_phase2','fn','time_grp','grp',...
    'wl','csc1','csc2','msgn','rrat','rdate')
keyboard

% -------------------------------------------------------------------------
function [R p] = dphist(time_grp1,rms_whisk,resp_phase2) 
edges = -pi:2*pi/18:pi;     % phase histogram bin limits
cnts = (edges(1:end-1) + edges(2:end)) / 2;     % phase histogram bin centers
sfwhisk = rms_whisk' .* time_grp1;
zeroinx = time_grp1==0;
spvr = [resp_phase2 sfwhisk];
spvr(zeroinx,:) = [];
spvr = sortrows(spvr,1);
[mn_phase mn_wh] = phasedep(spvr,edges);
figure
plot([cnts cnts+2*pi],[mn_wh mn_wh],'r')      % conditional distribution: E(whisk|phi1<phase<phi2)
if ~isempty(find(time_grp1,1))
    [R p] = lincirc_corr2(mn_wh',cnts')
else
    R = NaN;
    p = NaN;
end

% -------------------------------------------------------------------------
function [R p] = dphist2(time_grp1,time_grp2,rms_whisk,resp_phase2) 
edges = -pi:2*pi/18:pi;     % phase histogram bin limits
cnts = (edges(1:end-1) + edges(2:end)) / 2;     % phase histogram bin centers
sfwhisk = rms_whisk' .* time_grp1;
zeroinx = time_grp1==0;
spvr = [resp_phase2 sfwhisk];
spvr(zeroinx,:) = [];
spvr = sortrows(spvr,1);
[mn_phase mn_wh] = phasedep(spvr,edges);
figure
plot([cnts cnts+2*pi],[mn_wh mn_wh],'r')      % conditional distribution: E(whisk|phi1<phase<phi2)
if ~isempty(find(time_grp1,1))
    [R p] = lincirc_corr2(mn_wh',cnts')
else
    R = NaN;
    p = NaN;
end

sfwhisk = rms_whisk' .* time_grp2;
zeroinx = time_grp2==0;
spvr = [resp_phase2 sfwhisk];
spvr(zeroinx,:) = [];
spvr = sortrows(spvr,1);
[mn_phase mn_wh] = phasedep(spvr,edges);
hold on
plot([cnts cnts+2*pi],[mn_wh mn_wh],'k')     % overlay a second phase histogram