function afre_restrict_silicon_stand_phaseshift
%AFRE_RESTRICT_SILICON_STAND_PHASESHIFT   Phase analysis restricted to a frequency band.
%   AFRE_RESTRICT_SILICON_STAND_PHASESHIFT calculates and saves phase
%   analysis results for a given EEG frequency band (1.1 - 1.6 Hz) of
%   silicon probe data. EEG is standardized for phase calculations. See
%   APHASE_STAND for details on the analysis.
%
%   AFRE_RESTRICT_SILICON_STAND_PHASESHIFT performs a phase correction on
%   the EEG (see APHASESHIFT2 for details).
%
%   Reference:
%   Nelson MJ, Pouget P, Nilsen EA, Patten CD, Schall JD (2008) Review of
%   signal distortion through metal microelectrode recording circuits and 
%   filters. J Neurosci Methods 169:141-157. 
%
%   See also APHASESHIFT2 and AFRE_RESTRICT_SILICON_STAND.

% Input argument check
error(nargchk(0,0,nargin))

% Directories
global DATAPATH
inpdir = ['Y:\Names\Peter\balazs\sil2\'];
inpdir_phaseshift = [DATAPATH 'Andi\Ketxyl\Phaseshift\'];
resdir = [DATAPATH 'Andi\Ketxyl\FreBandRestrict_phase_silicon_stand_phaseshift\sil2\t2\'];
mm = pwd;

% Filelist
flist = {'d5.0' 'd5.2' 'd5.2b' 'd5.4' 'd5.5' 'd5.6' 'd5.6b' 'd5.7' 'd5.9'...
    'd6.0' 'd6.1' 'd6.2' 'd6.3' 'd6.4' 'd6.5' 'd6.6' 'd6.7'};     % for sil2
% flist = {'d5.5' 'd5.6' 'd5.7' 'd5.8' 'd6.0a' 'd6.0b' 'd6.2a' 'd6.2b'...
%     'd6.4a' 'd6.4b'};       % for sil3
sf = length(flist);

% Main
sr = 20900;
dsr = 1000;
const = sr / dsr;
edges = -180:20:180;     % edges for phase histogram
cnts = (edges(1:end-1) + edges(2:end)) / 2;
for o = 1:sf
    fname = flist{o}     % 'filename'
    ff = [inpdir fname '.mat'];       % load EEG
    load(ff)
    eeg = eeg1';
    len = length(eeg);
    clear eeg1 eeg2
    
    ff = [inpdir fname '.res.2'];       % load unit
    res = textread(ff);     % spike timestamps
    ff = [inpdir fname '.clu.2'];
    clu = textread(ff);     % indices of sorted cells
    ncells = clu(1) - 1;        % number of sorted cells + 1 (multiunit)
    clu = clu(2:end);
    
    dind = round(1:const:len);      % downsample EEG on 1000 Hz
    deeg = eeg(dind);
    
    [y,w] = b_fft3(deeg,1000);     % FFT
    inxs = w >= 0.21 & w <= 59.9;
    fnx = find(inxs);
    Hfft = figure;
    plot(w(inxs),abs(y(inxs)))
    load([inpdir_phaseshift 'phaseshift.mat'])  % phase correction
    lfnx = length(fnx);
    phsc = zeros(lfnx,1);
    for k = 1:lfnx
        phsc(k) = linterp(frs,phs,w(fnx(k)));
    end
    ng = angle(y);
    ng(fnx) = ng(fnx) + phsc / 180 * pi;
%     ng(fnx) = ng(fnx) + (phsc - 100) / 180 * pi;
    cf = abs(y) .* exp(i*ng);
    leny = length(y);
    if isequal(mod(leny,2),0)
        hlf = leny / 2;
        cf(hlf+2:end)=flipud(conj(cf(2:hlf)));
    else
        hlf = (leny + 1) / 2;
        cf(hlf+1:end)=flipud(conj(cf(2:hlf)));
    end
    ceeg = ifft(cf);
    deeg = real(ceeg)';
    fns = [fname '_FFT.fig'];     % save FFT figure
    cd(resdir)
    saveas(Hfft,fns)
    close(Hfft)
    clear Hfft
    
    nqf = dsr / 2;      % filtering EEG
    flt = fir1(4096,5/nqf,'low');      % lowpass filtering on 5 Hz
    feeg = filtfilt(flt,1,deeg);
    feeg = (feeg - mean(feeg)) / std(feeg);
    ahee = angle(hilbert(feeg));    % Hilbert-transformation
    
    for nc = 2:ncells
        vdisc = res(clu==nc)';
        seglen = 30 * sr;        % 30 sec. long segments
        lenr = floor(len/seglen);       % preallocation
        ind1 = [1:seglen:len];
        ind2 = ind1 + seglen -1;
        aang_as = [];
        for k = 1:lenr
            vd = vdisc(vdisc>ind1(k)&vdisc<ind2(k)) - ind1(k);      % localize
            lfeeg = feeg(round((ind1(k)-1)/const)+1:round(ind2(k)/const));
            lahee = ahee(round((ind1(k)-1)/const)+1:round(ind2(k)/const));

            vdisc2 = round(vd/const);
            
            cyclen = eegfre(lfeeg,lahee,dsr);    % slow cycle length
            freq = 1 / cyclen * 1000;
            if 1.1 < freq & freq < 1.6        % restrict frequency band
                [paang_as pdinx_as] = laphase_stand(lfeeg,lahee,vdisc2,dsr);    % PHASE - all spikes
                aang_as = [aang_as paang_as];
            end
        end
        if isempty(aang_as)     % continue if the freq. band is empty
            disp([fname ' The frequency band is empty.'])
            continue
        end
        H = figure;     % phase histogram
        n_as = length(aang_as);     % all spikes
        ftm_as = sum(exp(1).^(i*aang_as)) / n_as;    % first trigonometric moment
        ang_as = angle(ftm_as);   % mean angle
        mvl_as = abs(ftm_as);     % mean resultant length
        aang_as = aang_as * 180 / pi;
        ang_as = ang_as * 180 / pi;
        [nm_as,xout_as] = histc(aang_as,edges);   % phase histogram
        nm_as = nm_as(1:end-1);
        bar(cnts,nm_as')
        cmps = strread(fname,'%s','delimiter','_');
        titlestr = [];
        for tt = 1:length(cmps)
            titlestr = [titlestr ' ' cmps{tt}];
        end
        title(gca,[titlestr ' ' num2str(nc) ' all spikes'])
        x_lim = xlim;
        y_lim = ylim;
        axis([-200 200 y_lim(1) y_lim(2)])
        str = ['\it{Mean angle: }' '\bf ' num2str(ang_as)];
        text(-160,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','red')
        str = ['\it{Mean resultant length: }' '\bf ' num2str(mvl_as)];
        text(-160,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/12,str,'Color','red')
        str = ['\it{n: }' '\bf ' num2str(n_as)];
        text(-160,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/6,str,'Color','red')
        
        cd(resdir)
        fns = [fname '_' num2str(nc) '_PHASE.fig'];     % save figure
        saveas(H,fns)
        fn = [fname '_' num2str(nc) '_PHASE.mat'];      % save mat file
        save(fn,'aang_as','mvl_as');
        close(H)
        clear H
    end
end
cd(mm)

% -------------------------------------------------------------------------
function cyclen = eegfre(feeg,ahee,sr)

% Check SWS criteria:
% 1. discard cicles with EEG amp. lower then 2SD
% 2. discard cicles shorter then 100 ms (half wavelength of filter cutoff freq.)
fn = find(-diff(ahee)>2*pi-0.1);
sd = std(feeg);
cl6 = [];
for k = 1:length(fn)-1
    seeg = feeg(fn(k):fn(k+1));
    axs = max(seeg) - min(seeg);
    lg = (axs < 2 * sd)  | (fn(k+1) - fn(k) < 0.25 * sr);
    if ~lg
        cl6(end+1) = (fn(k+1) - fn(k)) / sr * 1000;   % remaining cycles' length in ms;
    end
end
cyclen = mean(cl6) / sr * 1000;   % cycle length in ms

% -------------------------------------------------------------------------
function [ang inx] = laphase_stand(feeg,ahee,vdisc,sr)
%LAPHASE_STAND    Phase angles for unit relative to EEG.
%   [A I] = LAPHASE_STAND(FEEG,AHEE,VDISC,SR) calculates Hilbert phase 
%   angles (A) for discriminated unit (VDISC) relative to filtered EEG 
%   (FEEG), when sampling frequency is given in SR and Hilbert-transform of
%   the EEG in AHEE. Cycles not fulfilling the following 2 criteria are
%   discarded: (i) EEG amp. higher then 2SD; (ii) min. 250 ms length.
%   Indices of discarded spikes of vdisc are returned in I.
%
%   See also HILBERT.

% Check SWS criteria:
% 1. discard cicles with EEG amp. lower then 2SD
% 2. discard cicles shorter then 250 ms
fn = find(-diff(ahee)>2*pi-0.1);
sd = std(feeg);
inx = find(vdisc<fn(1));
for k = 1:length(fn)-1
    seeg = feeg(fn(k):fn(k+1));
    axs = max(seeg) - min(seeg);
    sahee = ahee(fn(k):fn(k+1));
    if (axs < 2 * sd)  | (fn(k+1) - fn(k) < 0.25 * sr)
        inx = [inx find(vdisc>fn(k)&vdisc<fn(k+1))];
    end
end
inx = [inx find(vdisc>fn(end))];
vdisc(inx) = [];
ang = ahee(vdisc);