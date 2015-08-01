function afre_restrict_stand_layer5_w2w(inpdir1)
%AFRE_RESTRICT_STAND_LAYER5_W2W   Phase analysis restricted to a frequency band.
%   AFRE_RESTRICT_STAND_LAYER5_W2W(DR) calculates and saves wave-to-wave
%   phase analysis results for a given EEG frequency band (1.1 - 1.6 Hz).
%   Input directory should be given as an argument (DR). EEG is
%   standardized for phase calculations.
%
%   See also ACLUSTERCUT and APHASERUN_BURST2_STAND.

% Input argument check
error(nargchk(1,1,nargin))

% Directories
global DATAPATH
inpdir_bas = inpdir1;
inpdir2 = [DATAPATH 'Andi\Hajni_layer_5\Cluster\mat2\'];   % burst analysis data
resdir1 = [DATAPATH 'Andi\Hajni_layer_5\FreBandRestrict_phase_stand_w2w_c05\'];
resdir2 = [];
mm = pwd;

% Filelist
[files1_bas files_short1_bas] = filelist(inpdir_bas);
[files2 files_short2] = filelist2(inpdir2);
files_short_bas = intersect(files_short1_bas,files_short2);
sf_bas = length(files_short_bas);

% Progress indicator
[wb,awb1,awb2] = waitbar2([0 0],'Running AFRE RESTRICT STAND W2W...');
global WB
WB(end+1) = wb;

% Main
main(inpdir1,inpdir2,resdir1,resdir2,files_short_bas,sf_bas,wb,'bas');

close(wb)
cd(mm)

% -------------------------------------------------------------------------
function main(inpdir1,inpdir2,resdir1,resdir2,files_short,sf,wb,bob)

sr = 20000;
dsr = 1000;
const = sr / dsr;
for o = 1:sf
    fname = files_short{o}     % filename
    cmps = strread(fname(1:end-4),'%s','delimiter','_');     % waitbar
    if length(cmps) < 3
        strw = [cmps{1} ' ' cmps{2}];
    else
        strw = [cmps{1} ' ' cmps{2} ' ' cmps{3}];
    end
    waitbar2([(o-1)/sf 0],wb,strw);
    ff = [inpdir1 fname];       % load
    load(ff)
    len = length(eeg);
    clear data eeg0
    ff2 = [inpdir2 fname(1:end-6) '_CLUST2.mat'];
    load(ff2)
    
    nqf = dsr / 2;      % filtering EEG
    flt = fir1(4096,5/nqf,'low');      % lowpass filtering on 5 Hz
    feeg = filtfilt(flt,1,eeg(1:const:end));
    feeg = (feeg - mean(feeg)) / std(feeg);
    ahee = angle(hilbert(feeg));    % Hilbert-transformation
    
    fst = Burst(1,:);   % burst first spikes
    vb1 = vdisc(fst);
    sso = vdisc;       % single spikes
    ssi = vdisc;       % allfirstspikes
    ib = zeros(size(vdisc));       % intraburst spikes
    for k = 1:size(Burst,2)
        sso(Burst(1,k):Burst(2,k)) = 0;
        ssi(Burst(1,k)+1:Burst(2,k)) = 0;
        ib(Burst(1,k):Burst(2,k)) = 1;
    end
    sspo = sso(sso>0);
    afsp = ssi(ssi>0);
    ibs = vdisc(find(ib));
    vburst = vdisc(Burst);
    
    seglen = 30 * sr;        % 30 sec. long segments
    lenr = floor(len/seglen);       % preallocation
    ind1 = [1:seglen:len];
    ind2 = ind1 + seglen -1;
    aang_fs = [];
    aang_sp = [];
    aang_as = [];
    aang_afsp = [];
    ctafsp = [];
    aang_ibsang = [];
    aang_sspoang = [];
    aang_allang = [];
    ctallang = [];
    cycnb = [];
    r_fs = [];
    r_sp = [];
    r_as = [];
    r_afsp = [];
    H1 = figure;    % figures
    H2 = figure;
    for k = 1:lenr
        vd = vdisc(vdisc>ind1(k)&vdisc<ind2(k)) - ind1(k);      % localize
        loceeg = eeg(ind1(k):ind2(k));
        lfeeg = feeg((ind1(k)-1)/const+1:ind2(k)/const);
        lahee = ahee((ind1(k)-1)/const+1:ind2(k)/const);
        
% Phase histograms
        eeg2 = loceeg(1:const:end);    % downsample on 1000 Hz
        vdisc2 = round(vd/const);
        
        cyclen = eegfre(lfeeg,lahee,dsr);    % filter EEG, Hilbert-transform
        freq = 1 / cyclen * 1000;
        frlim1 = 1.1
        frlim2 = 1.6
        if frlim1 < freq & freq < frlim2        % restrict frequency band
            lvb1 = vb1(vb1>ind1(k)&vb1<ind2(k)) - ind1(k);
            lvb1 = round(lvb1/const);    % burst first spikes, downsample unit on 1000 Hz
%             [paang_fs pdinx_fs] = laphase_stand(lfeeg,lahee,lvb1,dsr);    % PHASE - burst first spikes
%             aang_fs = [aang_fs paang_fs];
%             ftm_fs0 = sum(exp(1).^(i*paang_fs)) / length(paang_fs);    % first trigonometric moment
%             mvl_fs0 = abs(ftm_fs0);     % mean resultant length
%             r_fs = [r_fs mvl_fs0];
%             
            lsspo = sspo(sspo>ind1(k)&sspo<ind2(k)) - ind1(k);
            lsspo = round(lsspo/const);    % single spikes, downsample unit on 1000 Hz
%             [paang_sp pdinx_sp] = laphase_stand(lfeeg,lahee,lsspo,dsr);    % PHASE - single spikes
%             aang_sp = [aang_sp paang_sp];
%             ftm_sp0 = sum(exp(1).^(i*paang_sp)) / length(paang_sp);    % first trigonometric moment
%             mvl_sp0 = abs(ftm_sp0);     % mean resultant length
%             r_sp = [r_sp mvl_sp0];
%                         
%             [paang_as pdinx_as] = laphase_stand(lfeeg,lahee,vdisc2,dsr);    % PHASE - all spikes
%             aang_as = [aang_as paang_as];
%             ftm_as0 = sum(exp(1).^(i*paang_as)) / length(paang_as);    % first trigonometric moment
%             mvl_as0 = abs(ftm_as0);     % mean resultant length
%             r_as = [r_as mvl_as0];
            
            lafsp = afsp(afsp>ind1(k)&afsp<ind2(k)) - ind1(k);
            lafsp = round(lafsp/const);    % all first spikes, downsample unit on 1000 Hz
            [paang_afsp pdinx_afsp] = laphase_stand(lfeeg,lahee,lafsp,dsr);    % PHASE - all first spikes
            aang_afsp = [aang_afsp paang_afsp];
            tafsp = afsp(afsp>ind1(k)&afsp<ind2(k));
            tafsp(pdinx_afsp) = [];
            ctafsp = [ctafsp tafsp];
            ftm_afsp0 = sum(exp(1).^(i*paang_afsp)) / length(paang_afsp);    % first trigonometric moment
            mvl_afsp0 = abs(ftm_afsp0);     % mean resultant length
            r_afsp = [r_afsp mvl_afsp0];
            
            libs = ibs(ibs>ind1(k)&ibs<ind2(k)) - ind1(k);
            libs = round(libs/const);    % intraburst spikes, downsample unit on 1000 Hz
            [paang_ibsang paang_sspoang paang_allang pcycnb pallang] = laphaseb(lfeeg,lahee,lsspo,libs,lvb1,dsr);
            aang_ibsang = [aang_ibsang paang_ibsang];       % PHASE - "cycle first"
            aang_sspoang = [aang_sspoang paang_sspoang];
            aang_allang = [aang_allang paang_allang];
            cycnb = [cycnb pcycnb];
            tallang = pallang + ind1(k) / const;
            ctallang = [ctallang tallang];
        end
        waitbar2([(o-1)/sf k/lenr],wb,strw);
    end
    if isempty(aang_afsp)
        close all
        if exist('fname')
            disp([fname ' The frequency band is empty.'])
        else
            disp(['No ' bob 'file.'])
        end
        continue
    end
    
% Wavelet
    eeg2 = eeg(1:const*5:end);    % downsample on 200 Hz
    [pow,phase,f] = eegwavelet(eeg2,200);        % WAVELET
    
% Plot
    cmps = strread(fname,'%s','delimiter','_');
    titlestr = [];
    for tt = 1:length(cmps)
        titlestr = [titlestr ' ' cmps{tt}];
    end
    figure(H1)
    set(gcf,'Position',[101 51 1150 750])     % maximize figure window
    subplot(2,1,1)
    plot(ctafsp/const/1000,aang_afsp*180/pi,'.')
    ylim([-180 180])
    xd = get(findobj(allchild(gca),'Type','line'),'XData');
    set(gca,'XLim',[min(xd) max(xd)])
    subplot(2,1,2)
    imagesc(pow)
    b_rescaleaxis('Y',f)
    ach = allchild(H1);     % figure title
    ax = findobj(ach,'type','axes');
    title(ax(end),titlestr)
    
    figure(H2)
    set(gcf,'Position',[101 51 1150 750])     % maximize figure window
    subplot(2,1,1)
    plot(ctallang/1000,aang_allang*180/pi,'.')
    ylim([-180 180])
    xd = get(findobj(allchild(gca),'Type','line'),'XData');
    set(gca,'XLim',[min(xd) max(xd)])
    subplot(2,1,2)
    imagesc(pow)
    b_rescaleaxis('Y',f)
    ach = allchild(H2);     % figure title
    ax = findobj(ach,'type','axes');
    title(ax(end),titlestr)
    
% Save
    cd(resdir1)
    fns = [fname(1:end-4) '_W2W_AFSP.fig'];
    saveas(H1,fns)
    fns = [fname(1:end-4) '_W2W_AFSP.jpg'];
    saveas(H1,fns)
    fns = [fname(1:end-4) '_W2W_CFSP.fig'];
    saveas(H2,fns)
    fns = [fname(1:end-4) '_W2W_CFSP.jpg'];
    saveas(H2,fns)
    
    fn = [fname(1:end-4) '_W2W.mat'];
    save(fn,'ctafsp','aang_afsp','ctallang','aang_allang');
    close all
end



% -------------------------------------------------------------------------
function [files2 files2_short] = filelist(inpdir)

% List of filenames
files = dir(inpdir);
files = files(3:end);
vrs = version;
if isequal(vrs(1:5),'7.4.0')
    files2 = struct('name',[],'date',[],'bytes',[],'isdir',[],'datenum',[]);
else
    files2 = struct('name',[],'date',[],'bytes',[],'isdir',[]);
end
files2_short = {};
for i = 1:length(files)
    if ~files(i).isdir
        files2(end+1) = files(i);
        files2_short{end+1} = files(i).name;
    end
end
files2 = files2(2:end);

% -------------------------------------------------------------------------
function [files2 files2_short] = filelist2(inpdir)

% List of filenames
files = dir(inpdir);
files = files(3:end);
vrs = version;
if isequal(vrs(1:5),'7.4.0')
    files2 = struct('name',[],'date',[],'bytes',[],'isdir',[],'datenum',[]);
else
    files2 = struct('name',[],'date',[],'bytes',[],'isdir',[]);
end
files2_short = {};
for i = 1:length(files)
    if ~files(i).isdir
        files2(end+1) = files(i);
        files2_short{end+1} = [files(i).name(1:end-11) '_d.mat'];
    end
end
files2 = files2(2:end);

% -------------------------------------------------------------------------
function [lvburst lburst] = locvburst(vburst,Burst,ind1,ind2)

fst = vburst(1,:);
lst = vburst(2,:);
fi = find(fst>ind1,1,'first');
li = find(lst<ind2,1,'last');
lvburst = vburst(:,fi:li);
lburst = Burst(:,fi:li);

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

% -------------------------------------------------------------------------
function [ang_ibsang ang_sspoang ang_allang cycnb allang] = laphaseb(feeg,ahee,sspo,ibs,vb1,sr)
%APHASEB    Phase angles for unit relative to EEG.
%   [IBSANG SSPOANG] = APHASEB(FEEG,AHEE,SSPO,IBS,VB1,SR) calculates
%   Hilbert phase angles for first intraburst and first single spike of
%   each cycle relative to the filtered EEG (FEEG), when sampling frequency
%   is given in SR, single spikes are given in SSPO, intraburst spikes are
%   given in IBS, burst first spikes are given in VB1 and Hilbert-transform
%   of the EEG in AHEE. Cycles not fulfilling the following 2 criteria are
%   discarded: (i) EEG amp. higher then 2SD; (ii) min. 250 ms length. 
%
%   [IBSANG SSPOANG ALLANG] = APHASEB(FEEG,AHEE,SSPO,IBS,VB1,SR) returns
%   the phase of cycle first spikes as well.
%
%   [IBSANG SSPOANG ALLANG CYCNB] = APHASEB(FEEG,AHEE,SSPO,IBS,VB1,SR)
%   returns the number of bursts in each cycles.
%
%   See also HILBERT.

% Check SWS criteria:
% 1. discard cicles with EEG amp. lower then 2SD
% 2. discard cicles shorter then 100 ms (half wavelength of filter cutoff freq.)
fn = find(-diff(ahee)>2*pi-0.1);
cyclen1 = mean(diff(fn)) / sr * 1000;   % cycle length in ms
sd = std(feeg);
ibsang = [];
sspoang = [];
allang = [];
cycnb = [];
for k = 1:length(fn)-1
    seeg = feeg(fn(k):fn(k+1));
    axs = max(seeg) - min(seeg);
    if ~(axs < 2 * sd)  | (fn(k+1) - fn(k) < 0.25 * sr)
        seeg = feeg(fn(k):fn(k+1));
        axs = max(seeg) - min(seeg);
        sahee = ahee(fn(k):fn(k+1));
        cycibs = ibs(ibs>fn(k)&ibs<fn(k+1));
        cycsspo = sspo(sspo>fn(k)&sspo<fn(k+1));
        cycvb1 = vb1(vb1>fn(k)&vb1<fn(k+1));
        if ~isempty(cycibs)
            if ~isempty(cycsspo)
                if cycibs(1) < cycsspo(1)
                    ibsang = [ibsang cycibs(1)];
                end
            else
                ibsang = [ibsang cycibs(1)];
            end
        end
        if ~isempty(cycsspo)
            if ~isempty(cycibs)
                if cycsspo(1) < cycibs(1)
                    sspoang = [sspoang cycsspo(1)];
                end
            else
                sspoang = [sspoang cycsspo(1)];
            end
        end
        if ~(isempty(cycibs) || isempty(cycsspo))
            allang = [allang min(cycibs(1),cycsspo(1))];
        elseif ~isempty(cycibs) && isempty(cycsspo)
            allang = [allang min(cycibs(1))];
        elseif isempty(cycibs) && ~isempty(cycsspo)
            allang = [allang min(cycsspo(1))];
        end
        cycnb(end+1) = length(cycvb1);
    end
end
ang_ibsang = ahee(ibsang);
ang_sspoang = ahee(sspoang);
ang_allang = ahee(allang);

% -------------------------------------------------------------------------
function st = asta(vdisc,eeg,wn)
%ASTA    Spike Triggered Average.
%   S = ASTA(VD,EEG,WN) calculates spike triggered EEG sum of EEG and
%   discriminated unit (VD) using WN windowsize. Each EEG window is
%   standardized.

% Standardize EEG
eeg = (eeg - mean(eeg)) / std(eeg);

% Calculate STA
wn2 = round(wn/2);
vdisc = vdisc(find(vdisc-wn2>0&vdisc+wn2<=length(eeg)));
lenv = length(vdisc);
st = zeros(lenv,2*wn2+1);
for t = 1:lenv
    if vdisc(t) - wn2 >= 1 & vdisc(t) + wn2 <= length(eeg)
        eeg2 = eeg(vdisc(t)-wn2:vdisc(t)+wn2);
        st(t,1:2*wn2+1) = eeg2;
    end
end