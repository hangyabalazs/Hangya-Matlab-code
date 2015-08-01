function lfwavelet
%LFWAVELET   Runs wavelet on a sequence of files.
%   LFWAVELET calculates and saves wavelet.
%
%   See also .

% Input argument check
error(nargchk(0,0,nargin))

% Directories
global DATAPATH
global DATADIR
inpdir = [DATADIR 'LentiFenti\'];
resdir = [DATAPATH 'LentiFenti\'];
mm = pwd;

% Filelist
[files files_short] = filelist(inpdir);
sf = length(files_short);

% Progress indicator
wb = waitbar(0,'Running LFWAVELET...','Position',[360 250 275 50]);
global WB
WB(end+1) = wb;

% Main
main(inpdir,resdir,files_short,sf,wb);

close(wb)
cd(mm)

% -------------------------------------------------------------------------
function main(inpdir,resdir,files_short,sf,wb)

sr = 20000;
dsr = 200;
const = sr / dsr;
nqf = dsr / 2;
for o = 1:sf
    fname = files_short{o}     % filename
    ff = [inpdir fname];       % load
    load(ff)
    eeg = MP_IL_S1.values';
    unit = Unit_row.values;
    
    eeg = eeg(1:const:end);     % downsampling EEG
    len = length(eeg);
    [im_eeg,f] = eegwavelet(eeg,dsr);        % EEG wavelet
    dt = 1 / dsr;
    wavetime = (0:len-1) * dt;
    
    H = figure;         % plot
    imagesc(im_eeg)
    mx = max(im_eeg(:));
    x_lim = get(gca,'XLim');
    y_lim = get(gca,'YLim');
    xcoord = 3 * ((x_lim(2)-x_lim(1)) + x_lim(1)) / 4;
    ycoord = 1 * ((y_lim(2)-y_lim(1)) + y_lim(1)) / 4;
    text(xcoord,ycoord,num2str(mx),'Color',[1 1 1]);
    ff = round(f*100) / 100;
    time = round(wavetime*100) / 100;
    b_rescaleaxis('Y',ff)
    b_rescaleaxis('X',time)
    tt = fname(1:end-4);
    tt(tt=='_') = ' ';
    tt = [tt ' EEG wavelet'];
    title(tt)
    setappdata(gca,'scalex',wavetime)
    setappdata(gca,'scaley',f)
    b_zoomset_for_wavelet
    fn = [resdir tt '.tiff'];       % save
%     saveas(H,fn)
    
    % Filter EEG
    flt = fir1(4096*4,1.5/nqf,'low');      % lowpass filtering on 5 Hz
    feeg = filtfilt(flt,1,eeg);
    feeg = (feeg - mean(feeg)) / std(feeg);
    ahee = angle(hilbert(feeg));    % Hilbert-transformation
    
    % Discrimination
    vdisc = disc(unit);
    vdisc = round(vdisc/const);   % downsample
    
    % Phase analysis
    seglen = 30 * dsr;        % 30 sec. long segments
    lenr = floor(len/seglen);       % preallocation
    ind1 = 1:seglen:len;
    ind2 = ind1 + seglen -1;
    for k = 1:lenr
        lvdisc = vdisc(vdisc>ind1(k)&vdisc<ind2(k)) - ind1(k);      % localize
        loceeg = eeg(ind1(k):ind2(k));
        lfeeg = feeg((ind1(k)-1)/const+1:ind2(k)/const);
        lahee = ahee((ind1(k)-1)/const+1:ind2(k)/const);
        
        % Phase histograms
        cyclen = eegfre(lfeeg,lahee,dsr);    % filter EEG, Hilbert-transform
        freq = 1 / cyclen * 1000
        if 2 < freq & freq < 2.5        % restrict frequency band
            llvdisc = lvdisc(lvdisc>ind1(k)&lvdisc<ind2(k)) - ind1(k);
            
            figure
            plot(loceeg,'g')
            hold on
            [paang_fs pdinx_fs] = laphase_stand(lfeeg,lahee,llvdisc,dsr);    % PHASE - burst first spikes
            aang_fs = [aang_fs paang_fs];
            ftm_fs0 = sum(exp(1).^(i*paang_fs)) / length(paang_fs);    % first trigonometric moment
            mvl_fs0 = abs(ftm_fs0);     % mean resultant length
            r_fs = [r_fs mvl_fs0];
            
            lsspo = sspo(sspo>ind1(k)&sspo<ind2(k)) - ind1(k);
            lsspo = round(lsspo/const);    % single spikes, downsample unit on 1000 Hz
            [paang_sp pdinx_sp] = laphase_stand(lfeeg,lahee,lsspo,dsr);    % PHASE - single spikes
            aang_sp = [aang_sp paang_sp];
            ftm_sp0 = sum(exp(1).^(i*paang_sp)) / length(paang_sp);    % first trigonometric moment
            mvl_sp0 = abs(ftm_sp0);     % mean resultant length
            r_sp = [r_sp mvl_sp0];
                        
            [paang_as pdinx_as] = laphase_stand(lfeeg,lahee,vdisc2,dsr);    % PHASE - all spikes
            aang_as = [aang_as paang_as];
            ftm_as0 = sum(exp(1).^(i*paang_as)) / length(paang_as);    % first trigonometric moment
            mvl_as0 = abs(ftm_as0);     % mean resultant length
            r_as = [r_as mvl_as0];
            
            lafsp = afsp(afsp>ind1(k)&afsp<ind2(k)) - ind1(k);
            lafsp = round(lafsp/const);    % all first spikes, downsample unit on 1000 Hz
            [paang_afsp pdinx_afsp] = laphase_stand(lfeeg,lahee,lafsp,dsr);    % PHASE - all first spikes
            aang_afsp = [aang_afsp paang_afsp];
            ftm_afsp0 = sum(exp(1).^(i*paang_afsp)) / length(paang_afsp);    % first trigonometric moment
            mvl_afsp0 = abs(ftm_afsp0);     % mean resultant length
            r_afsp = [r_afsp mvl_afsp0];
            
            libs = ibs(ibs>ind1(k)&ibs<ind2(k)) - ind1(k);
            libs = round(libs/const);    % intraburst spikes, downsample unit on 1000 Hz
            [paang_ibsang paang_sspoang paang_allang pcycnb] = laphaseb(lfeeg,lahee,lsspo,libs,lvb1,dsr);
            aang_ibsang = [aang_ibsang paang_ibsang];       % PHASE - "cycle first"
            aang_sspoang = [aang_sspoang paang_sspoang];
            aang_allang = [aang_allang paang_allang];
            cycnb = [cycnb pcycnb];
            
            wn = 2 * dsr;    % 2 sec. window
            st_fs = asta(lvb1,eeg2,wn);    % STA - burst first spikes
            st_sp = asta(lsspo,eeg2,wn);    % STA - single spikes
            st_as = asta(vdisc2,eeg2,wn);    % STA - all spikes
            st_afsp = asta(lafsp,eeg2,wn);    % STA - all first spikes
            nn_fs = nn_fs + length(lvb1);
            nn_sp = nn_sp + length(lsspo);
            nn_as = nn_as + length(vdisc2);
            nn_afsp = nn_afsp + length(lafsp);
        end
    end

    waitbar(o/sf)
    close all
end



% -------------------------------------------------------------------------
function [files2 files2_short] = filelist(inpdir)

% List of filenames
files = dir(inpdir);
files = files(3:end);
vrs = version;
if isequal(vrs(1),'7')
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
function [wave,f] = eegwavelet(dat,sr)

% Prepare for wavelet transformation
variance = std(dat) ^ 2;
dat = (dat - mean(dat)) / sqrt(variance) ;
n = length(dat);
dt = 1 / sr;
pad = 1;
dj = 0.02;    
j1 = ceil((1/dj) * log2(n/2));
j1 = ceil(j1);
j = (0:j1);
s0 = 2 * dt; 
s = s0 .* 2 .^ (j * dj);
omega0 = 6;
c = 4 * pi / (omega0 + sqrt(2+omega0^2));
fperiod = c .* s;
f = 1 ./ fperiod;
lag1 = 0.72;
param = -1;
mif = 0.5;          %minimal intresting frequency
mis = find(f>mif);
mis = mis(end);     %maximal intristing scale
mother = 'Morlet';

% Wavelet transformation
[wave,period,scale,coi] = b_wavelet_pow3(dat,dt,pad,dj,s0,j1,mother,param,mis);

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
plot(feeg)
plot(ahee,'k')
for k = 1:length(fn)-1
    seeg = feeg(fn(k):fn(k+1));
    axs = max(seeg) - min(seeg);
    sahee = ahee(fn(k):fn(k+1));
    if (axs < 2 * sd)  | (fn(k+1) - fn(k) < 0.25 * sr)
        inx = [inx find(vdisc>fn(k)&vdisc<fn(k+1))];
        plot(fn(k):fn(k+1),seeg,'r')
        plot(fn(k):fn(k+1),sahee,'r')
    end
end
inx = [inx find(vdisc>fn(end))];
vdisc(inx) = [];
ang = ahee(vdisc);