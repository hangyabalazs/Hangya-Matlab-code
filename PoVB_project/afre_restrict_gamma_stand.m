function afre_restrict_gamma_stand(inpdir1)
%AFRE_RESTRICT_GAMMA_STAND   Beta phase analysis restricted to a frequency band.
%   AFRE_RESTRICT_GAMMA_STAND(DR) calculates and saves beta phase (10 - 20
%   Hz) analysis results for a given EEG frequency band (1.1 - 1.6 Hz). See
%   APHASE_GAMMA for details on the analysis. Input directory should be
%   given as an argument (DR). EEG is standardized for phase calculations.
%
%   See also APHASE_GAMMA and AFRERESTRICTGAMMASTAND_CALL.

% Input argument check
error(nargchk(1,1,nargin))

% Directories
global DATAPATH
inpdir_bas = [inpdir1 'bas\'];
inpdir_bic = [inpdir1 'bic\'];
inpdir2 = [DATAPATH 'Andi\Ketxyl\Cluster\mat2\'];   % burst analysis data
resdir1 = [DATAPATH 'Andi\Ketxyl\FreBandRestrict_phase_spindle_stand\'];
mm = pwd;

% Filelist
[files1_bas files_short1_bas] = filelist(inpdir_bas);
[files1_bic files_short1_bic] = filelist(inpdir_bic);
[files2 files_short2] = filelist2(inpdir2);
files_short_bas = intersect(files_short1_bas,files_short2);
files_short_bic = intersect(files_short1_bic,files_short2);
sf_bas = length(files_short_bas);
sf_bic = length(files_short_bic);

% Progress indicator
[wb,awb1,awb2] = waitbar2([0 0],'Running AFRE RESTRICT GAMMA STAND...');
global WB
WB(end+1) = wb;

% Main
main(inpdir1,inpdir2,resdir1,files_short_bas,sf_bas,wb,'bas');
main(inpdir1,inpdir2,resdir1,files_short_bic,sf_bic,wb,'bic');

close(wb)
cd(mm)

% -------------------------------------------------------------------------
function main(inpdir1,inpdir2,resdir1,files_short,sf,wb,bob);

sr = 20000;
dsr = 1000;
const = sr / dsr;
edges = -180:20:180;     % edges for phase histogram
cnts = (edges(1:end-1) + edges(2:end)) / 2;
aang_fs = [];
aang_sp = [];
aang_as = [];
aang_afsp = [];
aang_ibsang = [];
aang_sspoang = [];
aang_allang = [];
cycnb = [];
betapower = [];
betapower_norm = [];
H3 = figure;
H5 = figure;
H6 = figure;
for o = 1:sf
    fname = files_short{o}     % filename
    cmps = strread(fname(1:end-4),'%s','delimiter','_');     % waitbar
    if length(cmps) < 3
        strw = [cmps{1} ' ' cmps{2}];
    else
        strw = [cmps{1} ' ' cmps{2} ' ' cmps{3}];
    end
    waitbar2([(o-1)/sf 0],wb,strw);
    ff = [inpdir1 bob '\' fname];       % load
    load(ff)
    eeg = data(:,2)';
    len = length(data);
    clear data eeg0
    ff2 = [inpdir2 fname(1:end-4) '_CLUST2.mat'];
    load(ff2)
    
    nqf = dsr / 2;      % filtering EEG
    flt = fir1(4096,5/nqf,'low');      % lowpass filtering on 5 Hz
    feeg_slow = filtfilt(flt,1,eeg(1:const:end));
    feeg_slow = (feeg_slow - mean(feeg_slow)) / std(feeg_slow);
    ahee_slow = angle(hilbert(feeg_slow));    % Hilbert-transformation
    flt = fir1(4096,[7 20]/nqf,'band');      % bandpass filtering
    feeg_beta = filtfilt(flt,1,eeg(1:const:end));
    feeg = (feeg_beta - mean(feeg_beta)) / std(feeg_beta);
    ahee_beta = angle(hilbert(feeg_beta));    % Hilbert-transformation
    
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
    for k = 1:lenr
        vd = vdisc(vdisc>ind1(k)&vdisc<ind2(k)) - ind1(k);      % localize
        loceeg = eeg(ind1(k):ind2(k));
        lfeeg_slow = feeg_slow((ind1(k)-1)/const+1:ind2(k)/const);
        lahee_slow = ahee_slow((ind1(k)-1)/const+1:ind2(k)/const);
        lfeeg_beta = feeg_beta((ind1(k)-1)/const+1:ind2(k)/const);
        lahee_beta = ahee_beta((ind1(k)-1)/const+1:ind2(k)/const);
        
% Phase histograms
        eeg2 = loceeg(1:const:end);    % downsample on 1000 Hz
        vdisc2 = round(vd/const);
        
        cyclen = eegfre(lfeeg_slow,lahee_slow,dsr);    % filter EEG, Hilbert-transform
        freq = 1 / cyclen * 1000;
        if 1.1 < freq & freq < 1.6        % restrict frequency band
            lvb1 = vb1(vb1>ind1(k)&vb1<ind2(k)) - ind1(k);
            lvb1 = round(lvb1/const);    % burst first spikes, downsample unit on 1000 Hz
            [paang_fs pdinx_fs] = laphase_stand(lfeeg_beta,lahee_beta,lvb1,dsr);    % PHASE - burst first spikes
            aang_fs = [aang_fs paang_fs];
            
            lsspo = sspo(sspo>ind1(k)&sspo<ind2(k)) - ind1(k);
            lsspo = round(lsspo/const);    % single spikes, downsample unit on 1000 Hz
            [paang_sp pdinx_sp] = laphase_stand(lfeeg_beta,lahee_beta,lsspo,dsr);    % PHASE - single spikes
            aang_sp = [aang_sp paang_sp];
                        
            [paang_as pdinx_as] = laphase_stand(lfeeg_beta,lahee_beta,vdisc2,dsr);    % PHASE - all spikes
            aang_as = [aang_as paang_as];
            
            lafsp = afsp(afsp>ind1(k)&afsp<ind2(k)) - ind1(k);
            lafsp = round(lafsp/const);    % all first spikes, downsample unit on 1000 Hz
            [paang_afsp pdinx_afsp] = laphase_stand(lfeeg_beta,lahee_beta,lafsp,dsr);    % PHASE - all first spikes
            aang_afsp = [aang_afsp paang_afsp];
            
            libs = ibs(ibs>ind1(k)&ibs<ind2(k)) - ind1(k);
            libs = round(libs/const);    % intraburst spikes, downsample unit on 1000 Hz
            [paang_ibsang paang_sspoang paang_allang pcycnb] = laphaseb(lfeeg_beta,lahee_beta,lsspo,libs,lvb1,dsr);
            aang_ibsang = [aang_ibsang paang_ibsang];       % PHASE - "cycle first"
            aang_sspoang = [aang_sspoang paang_sspoang];
            aang_allang = [aang_allang paang_allang];
            cycnb = [cycnb pcycnb];
            
            eeg2_200 = eeg2(1:5:end);        % downsample on 200 Hz
            [y,w] = b_fft2(eeg2_200,200);     % FFT
            inxs = w > 10 & w < 20;
            inxs2 = w > 5 & w < 25;
            integrated_power = sum(y(inxs));
            allpower = sum(y(2:end));
            betapower =  [betapower integrated_power];
            betapower_norm = [betapower_norm integrated_power/allpower];
        end
        waitbar2([(o-1)/sf k/lenr],wb,strw);
    end
end
if isempty(aang_fs) | isempty(aang_sp)
    close all
    if exist('fname')
        disp([fname ' The frequency band is empty.'])
    else
        disp(['No ' bob 'file.'])
    end
    return
end

n_fs = length(aang_fs);     % burst first spikes
ftm_fs = sum(exp(1).^(i*aang_fs)) / n_fs;    % first trigonometric moment
ang_fs = angle(ftm_fs);   % mean angle
mvl_fs = abs(ftm_fs);     % mean resultant length
aang_fs = aang_fs * 180 / pi;
ang_fs = ang_fs * 180 / pi;
[nm_fs,xout_fs] = histc(aang_fs,edges);   % phase histogram
nm_fs = nm_fs(1:end-1);

n_sp = length(aang_sp);     % single spikes
ftm_sp = sum(exp(1).^(i*aang_sp)) / n_sp;    % first trigonometric moment
ang_sp = angle(ftm_sp);   % mean angle
mvl_sp = abs(ftm_sp);     % mean resultant length
aang_sp = aang_sp * 180 / pi;
ang_sp = ang_sp * 180 / pi;
[nm_sp,xout_sp] = histc(aang_sp,edges);   % phase histogram
nm_sp = nm_sp(1:end-1);

n_as = length(aang_as);     % all spikes
ftm_as = sum(exp(1).^(i*aang_as)) / n_as;    % first trigonometric moment
ang_as = angle(ftm_as);   % mean angle
mvl_as = abs(ftm_as);     % mean resultant length
aang_as = aang_as * 180 / pi;
ang_as = ang_as * 180 / pi;
[nm_as,xout_as] = histc(aang_as,edges);   % phase histogram
nm_as = nm_as(1:end-1);

n_afsp = length(aang_afsp);     % all first spikes
ftm_afsp = sum(exp(1).^(i*aang_afsp)) / n_afsp;    % first trigonometric moment
ang_afsp = angle(ftm_afsp);   % mean angle
mvl_afsp = abs(ftm_afsp);     % mean resultant length
aang_afsp = aang_afsp * 180 / pi;
ang_afsp = ang_afsp * 180 / pi;
[nm_afsp,xout_afsp] = histc(aang_afsp,edges);   % phase histogram
nm_afsp = nm_afsp(1:end-1);

n_ibsang = length(aang_ibsang);     % "cycle first"
ftm_ibsang = sum(exp(1).^(i*aang_ibsang)) / n_ibsang;    % first trigonometric moment
ang_ibsang = angle(ftm_ibsang);   % mean angle
mvl_ibsang = abs(ftm_ibsang);     % mean resultant length
aang_ibsang = aang_ibsang * 180 / pi;
ang_ibsang = ang_ibsang * 180 / pi;
[nm_ibsang,xout_ibsang] = histc(aang_ibsang,edges);   % phase histogram
nm_ibsang = nm_ibsang(1:end-1);
n_sspoang = length(aang_sspoang);
ftm_sspoang = sum(exp(1).^(i*aang_sspoang)) / n_sspoang;    % first trigonometric moment
ang_sspoang = angle(ftm_sspoang);   % mean angle
mvl_sspoang = abs(ftm_sspoang);     % mean resultant length
aang_sspoang = aang_sspoang * 180 / pi;
ang_sspoang = ang_sspoang * 180 / pi;
if ~isempty(aang_sspoang)
    [nm_sspoang,xout_sspoang] = histc(aang_sspoang,edges);   % phase histogram
    nm_sspoang = nm_sspoang(1:end-1);
else
    nm_sspoang =zeros(1,length(edges)-1);
    xout_sspoang = xout_ibsang;
end
n_allang = length(aang_allang);
ftm_allang = sum(exp(1).^(i*aang_allang)) / n_allang;    % first trigonometric moment
ang_allang = angle(ftm_allang);   % mean angle
mvl_allang = abs(ftm_allang);     % mean resultant length
aang_allang = aang_allang * 180 / pi;
ang_allang = ang_allang * 180 / pi;
[nm_allang,xout_allang] = histc(aang_allang,edges);   % phase histogram
nm_allang = nm_allang(1:end-1);
[nm_cycnb xout_cycnb] = hist(cycnb(cycnb>0),(1:10));     % distr. of burst no./cycle

cmps = strread(fname,'%s','delimiter','_');
titlestr = [];
for tt = 1:length(cmps)
    titlestr = [titlestr ' ' cmps{tt}];
end

figure(H3)      % phase histograms
set(gcf,'Position',[1 31 1278 920])
subplot(2,2,1)      % all spikes
bar(cnts,nm_as'/length(aang_as))
title(gca,[titlestr ' all spikes'])
y_lim = ylim;
axis([-200 200 y_lim(1) y_lim(2)])
str = ['\it{Mean angle: }' '\bf ' num2str(ang_as)];
text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','red')
str = ['\it{Mean resultant length: }' '\bf ' num2str(mvl_as)];
text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/12,str,'Color','red')
str = ['\it{n: }' '\bf ' num2str(n_as)];
text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/6,str,'Color','red')
subplot(2,2,2)      % all first spikes
bar(cnts,nm_afsp'/length(aang_afsp))
title(gca,[titlestr ' all first spikes'])
y_lim = ylim;
axis([-200 200 y_lim(1) y_lim(2)])
str = ['\it{Mean angle: }' '\bf ' num2str(ang_afsp)];
text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','red')
str = ['\it{Mean resultant length: }' '\bf ' num2str(mvl_afsp)];
text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/12,str,'Color','red')
str = ['\it{n: }' '\bf ' num2str(n_afsp)];
text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/6,str,'Color','red')
subplot(2,2,3)      % burst first spikes
bar(cnts,nm_fs'/length(aang_fs))
title(gca,[titlestr ' burst first spikes'])
y_lim = ylim;
axis([-200 200 y_lim(1) y_lim(2)])
str = ['\it{Mean angle: }' '\bf ' num2str(ang_fs)];
text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','red')
str = ['\it{Mean resultant length: }' '\bf ' num2str(mvl_fs)];
text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/12,str,'Color','red')
str = ['\it{n: }' '\bf ' num2str(n_fs)];
text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/6,str,'Color','red')
subplot(2,2,4)      % single spikes
bar(cnts,nm_sp'/length(aang_sp))
title(gca,[titlestr ' single spikes'])
y_lim = ylim;
axis([-200 200 y_lim(1) y_lim(2)])
str = ['\it{Mean angle: }' '\bf ' num2str(ang_sp)];
text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','red')
str = ['\it{Mean resultant length: }' '\bf ' num2str(mvl_sp)];
text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/12,str,'Color','red')
str = ['\it{n: }' '\bf ' num2str(n_sp)];
text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/6,str,'Color','red')

figure(H5)
H5S1 = subplot(2,1,1);
[mu,kappa,Value,p,Rsquare] = b_watson2(aang_afsp/180*pi,1);     % fit of a von Mises distribution
H5S2 = subplot(2,1,2);
[param,err] = b_watsontwofit2(aang_afsp/180*pi,1);    % fit of the mixture of two von Mises distributions
figure(H6)
H6S1 = subplot(2,1,1);
[mu,kappa,Value,p,Rsquare] = b_watson2(aang_afsp/180*pi,1);     % fit of a von Mises distribution
H6S2 = subplot(2,1,2);
[param_,err_] = b_watsontwofit(aang_afsp/180*pi,1);    % fit of the mixture of two von Mises distributions

figure(H5)      % fit of a von Mises distribution
set(gcf,'Position',[1 31 1278 920])     % maximize figure window
axes(H5S1)
title(gca,[titlestr ' von Mises fit on all first spike distr.'])
y_lim = ylim;
str = ['\it{error: }' '\bf ' num2str(Rsquare)];
text(0.1,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','red')
str = ['\it{mu: }' '\bf ' num2str(mu*180/pi)];
text(0.1,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/12,str,'Color','red')
%     str = ['\it{kappa: }' '\bf ' num2str(kappa)];
%     text(2.1,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/12,str,'Color','red')
str = ['\it{R: }' '\bf ' num2str(b_A1(kappa))];
text(6.1,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/12,str,'Color','red')

axes(H5S2)     % fit of the mixture of two von Mises distributions
title(gca,[titlestr ' von Mises mixture fit on all first spike distr.'])
y_lim = ylim;
str = ['\it{error: }' '\bf ' num2str(err)];
text(0.1,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','red')
str = ['\it{mu: }' '\bf ' num2str(param(1)*180/pi) '   ' num2str(param(2)*180/pi)];
text(0.1,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/12,str,'Color','red')
%     str = ['\it{kappa: }' '\bf ' num2str(param(3)) '   ' num2str(param(4))];
%     text(3.1,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/12,str,'Color','red')
str = ['\it{R: }' '\bf ' num2str(b_A1(param(3))) '   ' num2str(b_A1(param(4)))];
text(6.1,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/12,str,'Color','red')
str = ['\it{p: }' '\bf ' num2str(param(5))];
text(0.1,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/6,str,'Color','red')

figure(H6)      % fit of a von Mises distribution
set(gcf,'Position',[1 31 1278 920])     % maximize figure window
axes(H6S1)
title(gca,[titlestr ' von Mises fit'])
y_lim = ylim;
str = ['\it{error: }' '\bf ' num2str(Rsquare)];
text(0.1,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','red')
str = ['\it{mu: }' '\bf ' num2str(mu*180/pi)];
text(0.1,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/12,str,'Color','red')
%     str = ['\it{kappa: }' '\bf ' num2str(kappa)];
%     text(2.1,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/12,str,'Color','red')
str = ['\it{R: }' '\bf ' num2str(b_A1(kappa))];
text(6.1,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/12,str,'Color','red')

axes(H6S2)     % fit of the mixture of two von Mises distributions
title(gca,[titlestr ' von Mises mixture fit'])
y_lim = ylim;
str = ['\it{error: }' '\bf ' num2str(err_)];
text(0.1,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','red')
str = ['\it{mu: }' '\bf ' num2str(param_(1)*180/pi) '   ' num2str(param_(2)*180/pi)];
text(0.1,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/12,str,'Color','red')
%     str = ['\it{kappa: }' '\bf ' num2str(param_(3)) '   ' num2str(param_(4))];
%     text(3.1,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/12,str,'Color','red')
str = ['\it{R: }' '\bf ' num2str(b_A1(param_(3))) '   ' num2str(b_A1(param_(4)))];
text(6.1,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/12,str,'Color','red')
str = ['\it{p: }' '\bf ' num2str(param_(5))];
text(0.1,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/6,str,'Color','red')

% Save
cd(resdir1)
fns = [fname(1:end-4) '_PHASEHIST.fig'];
saveas(H3,fns)
fns = [fname(1:end-4) '_PHASEHIST.jpg'];
saveas(H3,fns)
saveas(H5,fns)
fns = [fname(1:end-4) '_VMFIT1.jpg'];
saveas(H5,fns)
fns = [fname(1:end-4) '_VMFIT2.fig'];
saveas(H6,fns)
fns = [fname(1:end-4) '_VMFIT2.jpg'];
saveas(H6,fns)
fn = [fname(1:end-4) '_PHASE.mat'];
save(fn,'aang_fs','aang_sp','aang_as','aang_afsp','aang_ibsang',...
    'aang_sspoang','aang_allang','mvl_fs','mvl_sp','mvl_as','mvl_afsp',...
    'mvl_ibsang','mvl_sspoang','mvl_allang','cycnb','betapower','betapower_norm');
close all



% -------------------------------------------------------------------------
function [files2 files2_short] = filelist(inpdir)

% List of filenames
files = dir(inpdir);
files = files(3:end);
files2 = struct('name',[],'date',[],'bytes',[],'isdir',[]);
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
files2 = struct('name',[],'date',[],'bytes',[],'isdir',[]);
files2_short = {};
for i = 1:length(files)
    if ~files(i).isdir
        files2(end+1) = files(i);
        files2_short{end+1} = [files(i).name(1:end-11) '.mat'];
    end
end
files2 = files2(2:end);

% -------------------------------------------------------------------------
function cyclen = eegfre(feeg,ahee,sr)

% Check SWS criteria:
% 1. discard cicles with EEG amp. lower then 2SD
% 2. discard cicles shorter then 25 ms
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
%   discarded: (i) EEG amp. higher then 2SD; (ii) min. 50 ms length.
%   Indices of discarded spikes of vdisc are returned in I.
%
%   See also HILBERT.

% Check SWS criteria:
% 1. discard cicles with EEG amp. lower then 2SD
% 2. discard cicles shorter then 50 ms
fn = find(-diff(ahee)>2*pi-0.3);
sd = std(feeg);
inx = find(vdisc<fn(1));
for k = 1:length(fn)-1
    seeg = feeg(fn(k):fn(k+1));
    axs = max(seeg) - min(seeg);
    sahee = ahee(fn(k):fn(k+1));
    if (axs < 2 * sd)  | (fn(k+1) - fn(k) < 0.050 * sr)
        inx = [inx find(vdisc>fn(k)&vdisc<fn(k+1))];
    end
end
inx = [inx find(vdisc>fn(end))];
vdisc(inx) = [];
ang = ahee(vdisc);

% -------------------------------------------------------------------------
function [ang_ibsang ang_sspoang ang_allang cycnb] = laphaseb(feeg,ahee,sspo,ibs,vb1,sr)
%APHASEB    Phase angles for unit relative to EEG.
%   [IBSANG SSPOANG] = APHASEB(FEEG,AHEE,SSPO,IBS,VB1,SR) calculates
%   Hilbert phase angles for first intraburst and first single spike of
%   each cycle relative to the filtered EEG (FEEG), when sampling frequency
%   is given in SR, single spikes are given in SSPO, intraburst spikes are
%   given in IBS, burst first spikes are given in VB1 and Hilbert-transform
%   of the EEG in AHEE. Cycles not fulfilling the following 2 criteria are
%   discarded: (i) EEG amp. higher then 2SD; (ii) min. 50 ms length. 
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
% 2. discard cicles shorter then 50 ms
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
    if ~(axs < 2 * sd)  | (fn(k+1) - fn(k) < 0.050 * sr)
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