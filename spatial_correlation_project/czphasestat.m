function czphasestat
%CZPHASESTAT   Population phase histograms for hippocampal interneurons.
%   CZPHASESTAT plots and saves summary theta phase histograms for
%   hippocampal interneurons. Theta phase was calculated via
%   Hilbert-transform of the EEG filtered between 4 and 12 Hz (see
%   CZPHASE). Random samples (n = 5000) from individual phase distributions
%   are pooled to create population phase distributions for interneurons
%   forming negatively or positively correlated pairs. A control group
%   containing interneurons for which no correlated pyramidal neuron pairs
%   were found
%   is also calculated.
%
%   See also CZPHASE and CZPHASEPOLAR.

% Directories
global DATAPATH
pn = 'pos';
inpdir_eeg = [DATAPATH 'Czurko\czurko_EEG\'];
inpdir = [DATAPATH 'Czurko\czurko_EEG\phase\resubmission\pos\'];
resdir = inpdir;
mm = pwd;
cd(resdir)

% Load & create pooled sample
xlsname = [inpdir_eeg 'EEG3c.xls'];
[ntz mtz] = xlsread(xlsname,pn);
sf = size(mtz,1);   % number of pairs
dsc = 5000;
allangs = [];
for o = 1:sf
    disp(o)
    pfn = [mtz{o,1} ' ' mtz{o,2}];
    pfn(pfn=='_') = ' ';
    fn = [inpdir pfn '_ANGS.mat'];
    load(fn)        % load indiviual phase samples
    
    if p_rayleigh > 0.001   % skip non-phase locked cells
        continue
    end
    
    n_angs = length(angs);      % pooled sample
    rp = randperm(n_angs);
    try
        angs_ds = angs(rp(1:dsc));
        allangs = [allangs angs_ds];
    catch
    end
end

% Circular statistics
[ftm, hang, hmvl] = mvlmn(allangs,'rad');
hSE = circular_SE(allangs,'rad');
[hL hU] = circular_conf(allangs,'rad');
hSD = b_circular_std(allangs);

% Plot
edges = -180:20:180;     % edges for phase histogram
cnts = (edges(1:end-1) + edges(2:end)) / 2;
[nm,xout] = histc(allangs*180/pi,edges);   % phase histogram
nm = nm(1:end-1);
    
H = figure;
B = bar(cnts,nm'/length(allangs));
set(B,'FaceColor',[0.16 0.38 0.27])
y_lim = ylim;
axis([-200 200 y_lim(1) y_lim(2)])
str = ['\it{Mean angle: }' '\bf ' num2str(hang*180/pi)];
text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','red')
str = ['\it{Mean resultant length: }' '\bf ' num2str(hmvl)];
text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/12,str,'Color','red')
str = ['\it{n: }' '\bf ' num2str(length(allangs))];
text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/6,str,'Color','red')
pts = ['summary_phasehist_' pn];
ts = pts;
ts(ts=='_') = ' ';
title(ts)

% Save
fns = [pts '.fig'];
saveas(H,fns)
fns = [pts '.mat'];
save(fns,'allangs','hang','hmvl','ftm')