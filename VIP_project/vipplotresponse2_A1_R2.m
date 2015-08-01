function vipplotresponse2_A1_R2
%VIPPLOTRESPONSE2_A1_R2   Putative PV and SOM cells.
%   VIPPLOTRESPONSE2_A1_R2 compares fast spiking (FS) and slow firing (SS)
%   as well as narrow spiking (NS) and wide spiking (WS) 'inhibited'
%   neurons. (For details on the grouping of neurons based on the effect of
%   VIP stimulation, see VIPPLOTRESPONSE_A1_FCN.) Parameters of the
%   inhibitory effect are compared by Mann-Whitney U-test.
%
%   See also VIPPLOTRESPONSE2_A1_FCN and VIPISINFLUENCED_3B.

% Load PSTH variables
% load('C:\Balazs\_analysis\VIP\A1_psth_variables_revision1.mat')
load('C:\Balazs\_analysis\VIP\A1_psth_variables.mat')

% Groups
frlim = 1;   % lower firing rate limit for detecting inhibition
tagged = logical(vldty) & (isact==2);   % tagged cells
inx_act = logical(vldty) & (isact==1) & (isact~=2);  % activated cells
inx_inh = logical(vldty) & (isinh) & (baseline>frlim) & (isact~=2);   % inhibited cells; firing rate > 1Hz criterion
activated = find(inx_act&~inx_inh);  % indices of activated only cells
inhibited = find(inx_inh&~inx_act);  % indices of inhibited only cells
ai = find(inx_act&inx_inh);   % indices of cells with dual effect
inx = activation_peak(ai) < inhibition_peak(ai);
activated_inhibited = sort(ai(inx));  % activated, then inhibited
inhibited_activated = ai(~inx);   % inhibited, then activated
activated = [activated; activated_inhibited];   % activated and activated-inhibited
inhibited = [inhibited; inhibited_activated];   % inhibited and inhibited-activated
tagged = find(tagged);   % indices of tagged cells
[~, iainx] = intersect(inhibited,inhibited_activated);   % indices of inh-act within inh

% Relative FR change of inhibited cells
finh = minvalue(inhibited) ./ baseline(inhibited);   % rel. FR change for inhibited cells
r = rand(1,sum(finh==0)) / 700;
finh(finh==0) = r + 0.001;   % jitter cells that are inhibited to FR=0

% NS and WS cells
inhibited_NS = inhibited(spike_width(inhibited)<=275);   % NS: <=275 us
inhibited_WS = inhibited(spike_width(inhibited)>275);   % WS: >275 us

% Calculate medians
mNS_start = median(inhibition_start(inhibited_NS));   % inhibition start
mWS_start = median(inhibition_start(inhibited_WS));
mNS_peak = median(inhibition_peak(inhibited_NS));   % inhibition peak
mWS_peak = median(inhibition_peak(inhibited_WS));
mNS_end = median(inhibition_end(inhibited_NS));   % inhibition end
mWS_end = median(inhibition_end(inhibited_WS));
mNS_time = median(inhibition_time(inhibited_NS));   % inhibition duration
mWS_time = median(inhibition_time(inhibited_WS));

% Mann-Whitney U-test
[p H] = ranksum(inhibition_start(inhibited_NS),inhibition_start(inhibited_WS))
[p H] = ranksum(inhibition_peak(inhibited_NS),inhibition_peak(inhibited_WS))
[p H] = ranksum(inhibition_end(inhibited_NS),inhibition_end(inhibited_WS))
[p H] = ranksum(inhibition_time(inhibited_NS),inhibition_time(inhibited_WS))

% Box-whisker plot
boxstat(inhibition_end(inhibited_NS),inhibition_end(inhibited_WS),'NS','WS')
title('Inhibition end')
boxstat(inhibition_time(inhibited_NS),inhibition_time(inhibited_WS),'NS','WS')
title('Inhibition duration')

% FS ans SS cells
inhibited_FS = inhibited(baseline(inhibited)>5);
inhibited_SS = inhibited(baseline(inhibited)<=5);

% Calculate medians
mFS_start = median(inhibition_start(inhibited_FS));   % inhibition start
mSS_start = median(inhibition_start(inhibited_SS));
mFS_peak = median(inhibition_peak(inhibited_FS));   % inhibition peak
mSS_peak = median(inhibition_peak(inhibited_SS));
mFS_end = median(inhibition_end(inhibited_FS));   % inhibition end
mSS_end = median(inhibition_end(inhibited_SS));
mFS_time = median(inhibition_time(inhibited_FS));   % inhibition duration
mSS_time = median(inhibition_time(inhibited_SS));

% Mann-Whitney U-test
[p H] = ranksum(inhibition_start(inhibited_FS),inhibition_start(inhibited_SS))
[p H] = ranksum(inhibition_peak(inhibited_FS),inhibition_peak(inhibited_SS))
[p H] = ranksum(inhibition_end(inhibited_FS),inhibition_end(inhibited_SS))
[p H] = ranksum(inhibition_time(inhibited_FS),inhibition_time(inhibited_SS))

% Box-whisker plot
boxstat(inhibition_end(inhibited_FS),inhibition_end(inhibited_SS),'FS','SS')
title('Inhibition end')
boxstat(inhibition_time(inhibited_FS),inhibition_time(inhibited_SS),'FS','SS')
title('Inhibition duration')