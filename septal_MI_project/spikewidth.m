function spikewidth
%SPIKEWIDTH   Spike width and firing rate calculation.
%   SPIKEWIDTH calculates spike width from positive to negative peak, mean
%   firing rate and peak firing rate calculated in 5-second-long windows. 

% Input argument check
error(nargchk(0,0,nargin))

% Directories
global DATADIR
global DATAPATH
inpdir = [DATADIR 'mi_zshift_data\hippo_cells\'];
inpdir2 = [DATADIR 'mi_zshift_data\discriminated_hc\'];
resdir = [DATAPATH 'Burst\Spikewidth\'];
mm = pwd;
dr = dir(inpdir);
dr = dr(3:end);
sf = length(dr);

% Import
sr = 10000;     % sampling rate
mspw = zeros(1,sf);     % mean spike width
mfr = zeros(1,sf);      % mean firing rate
mpfr = zeros(1,sf);      % mean peak firing rate
for k = 1:sf
    ff = [inpdir dr(k).name];   % load original data
    load(ff)
    if exist('data')
        unit = data(:,2)';
    else
        unit = Unit.values';
    end
    len = length(unit);
    ff2 = [inpdir2 dr(k).name(1:end-4) '_d.mat'];     % load discriminated data
    try
        load(ff2)
    catch
        continue
    end
    
% Spike width
    spw = zeros(1,length(vdisc));
    for t = 1:length(vdisc)
        psp = vdisc(t);     % spike peak
        msp = unit(vdisc(t):vdisc(t)+3*sr/1000);
        nsp = psp + find(msp==min(msp)) / sr * 1000;  % negative peak
        if nsp == length(msp)
            error('No negative peak.')
        end
        spw(t) = nsp(1) - psp;
    end
    mspw(k) = mean(spw);
    
% Firing rate
    efflen = (vdisc(end) - vdisc(1)) / sr;
    mfr(k) = (length(vdisc) - 1) / efflen;
    seglen = 5 * sr;        % 5 sec. long segments
    lenr = floor(len/seglen);       % preallocation
    ind1 = [1:seglen:len];
    ind2 = ind1 + seglen - 1;
    ppfr = zeros(1,lenr);
    for sm = 1:lenr
        vd = vdisc(vdisc>ind1(sm)&vdisc<ind2(sm)) - ind1(sm);      % localize
        efflen = seglen / sr;
        ppfr(sm) = (length(vd) - 1) / efflen;
    end
    mpfr(k) = max(ppfr);
end
mean(mspw)
mean(mfr)
mean(mpfr)