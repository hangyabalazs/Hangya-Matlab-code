function [mn_phase mn_rt mn mvl inx R p H Hc] = phasedep2(eeg,edges,depy,sr,fr1,fr2,fo,mlt,ip)
%PHASEDEP2   Phase dependence.
%   [MP MR MN MVL INX R P H HC] = PHASEDEP2(EEG,EDGES,S,SR,FR1,FR2,FO,MLT,IP)
%   calculates phase dependence of the variable given in S relative to EEG
%   sampled at SR with respect to phase bins determined by EDGES. Phase is
%   calculated via Hilbert-transform of the EEG filtered between FR1 and
%   FR2 Hz and standardized. Filtering is performed by FIR filter applying
%   bidirectional zero-phase filtering with a filter order of FO (FO,
%   optional input parameter; default, 2048). Output arguments: 
%       MR: mean values conditioned on the phase intervals
%       MP: sample sizes in the phase bins
%       MN: weighted mean of angles
%       MVL: mean length of weighted sum vector
%       INX: indeces of EEG left out from the calculations (see criteria 
%           below)
%       R: linear-circular correlation coefficient (for bin centers and 
%           mean values)
%       P: corresponding significance value
%       H: handle for output plot
%       HC: handle for control phase histogram (see below)
%   Cycles with a length out of the frequency range or with an amplitude
%   lower than mean + MLT * SD are discarded (MLT, optional input
%   parameter; default value, 2). IP, optional input argument; if IP equals
%   to 1, control plots are displayed (default, 0). Control plots include
%   filtered trace and phase series as well as phase histogram for all EEG
%   phase values.
%
%   See also PHASEHIST and PHASEDEP.

% Input argument check
error(nargchk(6,9,nargin))
if nargin < 9
    ip = 0;
end
if nargin < 8
    mlt = 2;
end
if nargin < 7
    fo = 2048;
end

% Filtering EEG
nqf = sr / 2;
flt = fir1(fo,[fr1 fr2]/nqf,'band');      % bandpass filtering
feeg = filtfilt(flt,1,eeg);
feeg = standardize(feeg);

% Hilbert transformation of the EEG
ahee = angle(hilbert(feeg));
aheedeg = ahee * (180 / pi);
if ip
    figure
    plot(standardize(eeg))
    hold on
    plot(feeg,'g')
    plot(ahee,'c')
end

% Check criteria:
% 1. discard cicles with EEG amp. lower then MLT * SD
% 2. discard cicles shorter then 1/fr2 s
% 3. discard cicles longer then 1/fr1 s
fn0 = valuecrossing(1:length(ahee),ahee',0,'down');
fn = round(fn0);
sd = std(feeg);
inx = [];
for k = 1:length(fn)-1
    seeg = feeg(fn(k)+1:fn(k+1));
    axs = max(seeg) - min(seeg);
    if (axs < mlt * sd)  || (fn(k+1) - fn(k) > 1 / fr1 * sr) || (fn(k+1) - fn(k) < 1 / fr2 * sr)
        inx = [inx fn(k)+1:fn(k+1)];
        if ip
            plot(fn(k)+1:fn(k+1),seeg,'r')
        end
    end
end
ahee(inx) = [];
depy(inx) = [];

% Phase dependence
n = length(edges);
mn_phase = zeros(1,n-1);
mn_rt = zeros(1,n-1);
for k = 2:n
    inx = find(ahee>edges(k-1)&ahee<edges(k));
    mn_phase(k-1) = length(inx);
    mn_rt(k-1) = mean(depy(inx));
end

% Mean
pmn = depy' .* exp(i*ahee);
cmn = sum(pmn) / length(pmn);
mn = angle(cmn);
mvl = abs(cmn) / mean(depy);

% Linear-circular correlation
cnts = (edges(1:end-1) + edges(2:end)) / 2;     % phase histogram bin centers
[R p] = lincirc_corr2(mn_rt',cnts');

% Plot phase dependence
H = figure;
plot([cnts*180/pi cnts*180/pi+360],[mn_rt mn_rt],'r')      % conditional distribution
y_lim = ylim;
ylim([0 y_lim(2)])
str = ['\it{Mean angle: }' '\bf ' num2str(mn*180/pi)];
text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','black')
str = ['\it{Mean resultant length: }' '\bf ' num2str(mvl)];
text(60,y_lim(2)-1.5*(y_lim(2)-y_lim(1))/6,str,'Color','black')
str = ['\it{n: }' '\bf ' num2str(length(ahee))];
text(60,y_lim(2)-2*(y_lim(2)-y_lim(1))/6,str,'Color','black')
if p < 0.05
    clr = 'red';
else
    clr = 'black';
end
str = ['\it{circ-lin corr R = }' '\bf ' num2str(R)];
text(60,y_lim(2)-2.5*(y_lim(2)-y_lim(1))/6,str,'Color',clr)
str = ['\it{circ-lin corr p < }' '\bf ' num2str(p)];
text(60,y_lim(2)-3*(y_lim(2)-y_lim(1))/6,str,'Color',clr)

% Plot phase histogram of all EEG phase values
if ip || nargout > 8
    [nm,xout] = histc(ahee*180/pi,edges*180/pi);   % control phase histogram (all EEG phase values)
    nm = nm(1:end-1);
    Hc = figure;
    B = bar(cnts*180/pi,nm'/length(ahee));
    set(B,'FaceColor',[0.16 0.38 0.27])
    y_lim = ylim;
    axis([-200 200 y_lim(1) y_lim(2)])
end