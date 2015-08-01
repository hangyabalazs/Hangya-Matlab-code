function afre_restrict_stand(inpdir1)
%AFRE_RESTRICT_STAND   Phase and burst analysis restricted to a frequency band.
%   AFRE_RESTRICT_STAND(DR) calculates and saves phase and burst analysis
%   results for a given EEG frequency band (1.1 - 1.6 Hz). See ACLUSTERCUT
%   and APHASE_BURST2_STAND for details on the analysis. Input directory
%   should be given as an argument (DR). EEG is standardized for phase
%   calculations.
%
%   See also ACLUSTERCUT and APHASERUN_BURST2_STAND.

% Input argument check
error(nargchk(1,1,nargin))

% Directories
global DATAPATH
inpdir_bas = [inpdir1 'bas\'];
inpdir_bic = [inpdir1 'bic\'];
% inpdir2 = [DATAPATH 'Andi\Hajni_layer_5\Cluster\mat2\'];   % burst analysis data
% resdir1 = [DATAPATH 'Andi\Hajni_layer_5\FreBandRestrict_phase_stand\'];
% resdir2 = [DATAPATH 'Andi\Hajni_layer_5\FreBandRestrict_burst_stand\'];
inpdir2 = [DATAPATH 'Andi\Ketxyl\Cluster\mat2\'];   % burst analysis data
resdir1 = [DATAPATH 'Andi\Ketxyl\RawDataFigs\'];
resdir2 = [DATAPATH 'Andi\Ketxyl\RawDataFigs\'];
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
[wb,awb1,awb2] = waitbar2([0 0],'Running AFRE RESTRICT STAND...');
global WB
WB(end+1) = wb;

% Main
main(inpdir1,inpdir2,resdir1,resdir2,files_short_bas,sf_bas,wb,'bas');
main(inpdir1,inpdir2,resdir1,resdir2,files_short_bic,sf_bic,wb,'bic');

close(wb)
cd(mm)

% -------------------------------------------------------------------------
function main(inpdir1,inpdir2,resdir1,resdir2,files_short,sf,wb,bob)

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
r_fs = [];
r_sp = [];
r_as = [];
r_afsp = [];
preburst = [];
burstspikeno = 0;
allspikeno = 0;
intraburstfreq = [];
ibspno = [];
burstlength = [];
bfreqnumer = 0;
bfreqdenom = 0;
efflen = 0;
fratenumer = 0;
isimtx = {};
H1ibspno = [];
H1ibfr = [];
H1bl = [];
H1aang_fs = [];
nn_fs = 0;
nn_sp = 0;
nn_as = 0;
nn_afsp = 0;
H1 = figure;    % figures
H2 = figure;
H3 = figure;
H4 = figure;
H5 = figure;
H6 = figure;
H8 = figure;
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
    for k = 1:lenr
        vd = vdisc(vdisc>ind1(k)&vdisc<ind2(k)) - ind1(k);      % localize
        loceeg = eeg(ind1(k):ind2(k));
        lfeeg = feeg((ind1(k)-1)/const+1:ind2(k)/const);
        lahee = ahee((ind1(k)-1)/const+1:ind2(k)/const);
        
% Phase histograms
        eeg2 = loceeg(1:const:end);    % downsample on 1000 Hz
        vdisc2 = round(vd/const);
        
        cyclen = eegfre(lfeeg,lahee,dsr);    % filter EEG, Hilbert-transform
        freq = 1 / cyclen * 1000
        if 2 < freq & freq < 2.5        % restrict frequency band
            lvb1 = vb1(vb1>ind1(k)&vb1<ind2(k)) - ind1(k);
            lvb1 = round(lvb1/const);    % burst first spikes, downsample unit on 1000 Hz
            figure
            plot(eeg2,'g')
            hold on
            [paang_fs pdinx_fs] = laphase_stand(lfeeg,lahee,lvb1,dsr);    % PHASE - burst first spikes
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
            
% Burst statistics
            [lvb lburst] = locvburst(vburst,Burst,ind1(k),ind2(k));
            preburst = [preburst lburst];
            burstnum = size(lvb,2);
            intraburstiv = [];
            intraburstnum = zeros(1,burstnum);
            for j = 1:burstnum      % computing intraburstiv
                b = vdisc(vdisc>=lvb(1,j)&vdisc<=lvb(2,j));
                db = diff(b);
                intraburstiv = [intraburstiv db];
                intraburstnum(j) = length(b);   %intraburst spike number
                for m = 1:length(b) - 1     % ISI matrix
                    if size(isimtx,1) < (length(b) - 1) | size(isimtx,2) < m
                        isimtx{length(b)-1,m} = [];
                    end
                    isimtx{length(b)-1,m} = [isimtx{length(b)-1,m} (vdisc(lburst(1,j)+m)-vdisc(lburst(1,j)+m-1))/20]; %in ms
                end
            end
            burstspikeno = burstspikeno + (length(intraburstiv) + burstnum);
            allspikeno = allspikeno + length(vd);
            burstlen = (lvb(2,:) - lvb(1,:)) / sr;
            burstlength = [burstlength burstlen];
            if ~isempty(intraburstnum)
                intraburstfreq = [intraburstfreq (intraburstnum-1)./burstlen];
            end
            ibspno = [ibspno intraburstnum];
            bfreqnumer = bfreqnumer + sr * (burstnum - 1);
            bfreqdenom = bfreqdenom + (lvb(2,end) - lvb(1,1));
            efflen = efflen + (vd(end) - vd(1)) / sr;
            fratenumer = fratenumer + length(vd) - 1;
            
            pH1ibspno = intraburstnum;     % for H1 fig.
            pH1ibfr = (intraburstnum-1)./burstlen;
            pH1bl = burstlen;
            if ~isempty(pdinx_fs) && pdinx_fs(end) > burstnum + 1
                error('Technical error 208.')
            end
            pH1ibspno(pdinx_fs(pdinx_fs<=burstnum)) = [];
            pH1ibfr(pdinx_fs(pdinx_fs<=burstnum)) = [];
            pH1bl(pdinx_fs(pdinx_fs<=burstnum)) = [];
            H1ibspno = [H1ibspno pH1ibspno];
            H1ibfr = [H1ibfr pH1ibfr];
            H1bl = [H1bl pH1bl];
            if length(lvb1) == size(lvb,2)       % 'vb1' may contain the first spike of a truncated burst
                pH1aang_fs = paang_fs;
            elseif length(lvb1) == size(lvb,2) + 1
                if ~isempty(pdinx_fs) && pdinx_fs(end) == burstnum + 1
                    pH1aang_fs = paang_fs;      % the truncated burst may already been taken out from the phase analysis
                else
                    pH1aang_fs = paang_fs(1:end-1);
                end
            else
                error('Matrix dimension mismatch.')
            end
            H1aang_fs = [H1aang_fs pH1aang_fs*180/pi];
        end
        waitbar2([(o-1)/sf k/lenr],wb,strw);
    end
end
if isempty(aang_fs)
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

Burstiness = burstspikeno / allspikeno;     % burst parameters
IntraBurstFrequency.mean = mean(intraburstfreq);
IntraBurstFrequency.sd = std(intraburstfreq);
IntraBurstFrequency.all = intraburstfreq;
IntraBurstSpikeNumber.mean = mean(ibspno);
IntraBurstSpikeNumber.sd = std(ibspno);
IntraBurstSpikeNumber.all = ibspno;
BurstLength.mean = mean(burstlength);
BurstLength.sd = std(burstlength);
BurstLength.all = burstlength;
BurstFrequency = bfreqnumer / bfreqdenom;
FiringRate = fratenumer / efflen;
Burst = preburst;
for x = 1:size(isimtx,1)
    warning off
    for y = 1:size(isimtx,2)
        IsiMatrix.mean(x,y) = mean(isimtx{x,y});
        IsiMatrix.sd(x,y) = std(isimtx{x,y});
        IsiMatrix.num(x,y) = length(isimtx{x,y});
    end
    warning backtrace
end

sta_fs = mean(st_fs,1);       % STA
sta_sp = mean(st_sp,1);
sta_as = mean(st_as,1);
sta_afsp = mean(st_afsp,1);
sta_index1_fs = max(sta_fs) - mean(sta_fs);
sta_index1_sp = max(sta_sp) - mean(sta_sp);
sta_index1_as = max(sta_as) - mean(sta_as);
sta_index1_afsp = max(sta_afsp) - mean(sta_afsp);
sta_index2_fs = max(sta_fs);
sta_index2_sp = max(sta_sp);
sta_index2_as = max(sta_as);
sta_index2_afsp = max(sta_afsp);

figure(H1)      % burst parameters vs. phase
set(gcf,'Position',[1 31 1278 920])     % maximize figure window
subplot(2,2,1)
% bar(xout,nm_fs/length(aang_fs))
bar(cnts,nm_afsp')
cmps = strread(fname,'%s','delimiter','_');
titlestr = [];
for tt = 1:length(cmps)
    titlestr = [titlestr ' ' cmps{tt}];
end
title(gca,[titlestr ' all first spikes'])
x_lim = xlim;
y_lim = ylim;
axis([-200 200 y_lim(1) y_lim(2)])
str = ['\it{Mean angle: }' '\bf ' num2str(ang_afsp)];
text(-160,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','red')
str = ['\it{Mean resultant length: }' '\bf ' num2str(mvl_afsp)];
text(-160,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/12,str,'Color','red')
str = ['\it{n: }' '\bf ' num2str(n_afsp)];
text(-160,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/6,str,'Color','red')

subplot(2,2,2)
plot(H1ibspno+rand(size(H1ibspno))*0.1,H1aang_fs,'.')     % intraburst spike number vs. phase
hold on
plot(ones(size(aang_sp))+rand(size(aang_sp))*1-1,aang_sp,'.')
title(gca,'intraburst spike number')
subplot(2,2,3)
plot(H1ibfr,H1aang_fs,'.','Color','white')
hold on
title(gca,'intraburst frequency')     % intraburst freq. vs. phase
clr = {'red' 'green' 'blue' 'magenta' 'yellow' 'cyan' ...
    [0 153 0]/256 [255 102 0]/256 [255 204 204]/256 [102 102 102]/256};
mibs = max(H1ibspno);
plot_handle = [];
legstr = {};
for k = 2:mibs
    str = ['let = find(H1ibspno==' num2str(k) ');'];
    eval(str)
    if ~isempty(let)
        plot_handle(end+1) = plot(H1ibfr(let),H1aang_fs(let),'.','Color',clr{k-1});
        legstr{end+1} = num2str(k);
    end
end
set(gca,'Color','black')
[L,OBJ] = legend(plot_handle,legstr,'Location','NorthWestOutside');
fo = findobj(OBJ,'Type','text');
set(fo,'Color','white')
subplot(2,2,4)
plot(H1bl,H1aang_fs,'.','Color','white')
hold on
title(gca,'burst length')     % burst length vs. phase
plot_handle = [];
legstr = {};
for k = 2:mibs
    str = ['let = find(H1ibspno==' num2str(k) ');'];
    eval(str)
    if ~isempty(let)
        plot_handle(end+1) = plot(H1bl(let),H1aang_fs(let),'.','Color',clr{k-1});
        legstr{end+1} = num2str(k);
    end
end
set(gca,'Color','black')
[L,OBJ] = legend(plot_handle,legstr,'Location','NorthEastOutside');
fo = findobj(OBJ,'Type','text');
set(fo,'Color','white')

figure(H2)      % burst parameters vs. phase - colorcoded & 3D
nmm = [];
xoutt = [];
legstr_ = {};
ftm = zeros(1,mibs-1);
mn_rad = zeros(1,mibs-1);
mn = zeros(1,mibs-1);
mrl = zeros(1,mibs-1);
ftm_ = zeros(1,mibs-1);
Lconf = zeros(1,mibs-1);
Uconf = zeros(1,mibs-1);
for k = 2:mibs
    str = ['let = find(H1ibspno==' num2str(k) ');'];
    eval(str)
    [nm,xout] = hist(H1aang_fs(let),edges);
    nmm = [nmm; nm];
    xoutt = [xoutt; xout(:)'];
    lele = length(let);
    ng = H1aang_fs(let) / 180 * pi;
    ftm(k-1) = sum(exp(1).^(i*ng)) / lele;
    if lele > 10
        mn_rad(k-1) = angle(ftm(k-1));
        mn(k-1) = angle(ftm(k-1)) * 180 / pi;
        mrl(k-1) = abs(ftm(k-1));
        ftm_(k-1) = ftm(k-1);
        [Lconf(k-1), Uconf(k-1)] = wconf(ng,mn_rad(k-1),mrl(k-1));
        legstr_{end+1} = num2str(k);
    else
        mn_rad(k-1) = NaN;
        mn(k-1) = NaN;
        mrl(k-1) = NaN;
        ftm_(k-1) = NaN;
        Lconf(k-1) = NaN;
        Uconf(k-1) = NaN;
    end
end
S = subplot(2,2,1);
CMP = compass(real(exp(1).^(i*0.1)),imag(exp(1).^(i*0.1)));
set(CMP,'Color','black')
hold on
CMP = [];
for k = 2:mibs
    if ~isnan(ftm_(k-1))
        CMP(end+1) = compass(real(ftm_(k-1)),imag(ftm_(k-1)));
        set(CMP(end),'Color',clr{k-1},'LineWidth',2)
    end
end
fob = findobj(allchild(S),'Type','patch');
set(fob,'FaceColor','black')
[L,OBJ] = legend(CMP,legstr_,'Location','NorthEastOutside');
fo = findobj(OBJ,'Type','text');
set(fo,'Color','white')
set(L,'Color','black')
subplot(2,2,2)
pcolor(edges,[2:mibs+1],[nmm; zeros(1,size(nmm,2))])
shading flat
colorbar
title(titlestr)
subplot(2,2,3)
bar3(xout',(flipud(nmm))',0.35);
view(-60,20)
%     jcm = colormap(jet);
%     colormap(flipud(jet))
set(gca,'XTickLabel',[2:mibs])
subplot(2,2,4)
plot([2:mibs],mn)
for k = 2:mibs
    if ~isnan(mn(k-1))
        text(k,mn(k-1)+5,num2str(mrl(k-1)))
        line([k k],[Lconf(k-1)/pi*180 Uconf(k-1)/pi*180],'Color','red')
    end
end
xlim([1 mibs])
set(gcf,'Position',[1 31 1278 920])     % maximize figure window

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

figure(H8)      % "cycle first" phase histograms
set(gcf,'Position',[1 31 1278 920])
subplot(2,2,1)      % all cycle first spikes
bar(cnts,nm_allang'/length(aang_allang))
title(gca,[titlestr ' all cycle first spikes'])
y_lim = ylim;
axis([-200 200 y_lim(1) y_lim(2)])
str = ['\it{Mean angle: }' '\bf ' num2str(ang_allang)];
text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','red')
str = ['\it{Mean resultant length: }' '\bf ' num2str(mvl_allang)];
text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/12,str,'Color','red')
str = ['\it{n: }' '\bf ' num2str(n_allang)];
text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/6,str,'Color','red')
subplot(2,2,2)      % cycle first intraburst spikes
bar(cnts,nm_ibsang'/length(aang_ibsang))
title(gca,[titlestr ' cycle first intraburst spikes'])
y_lim = ylim;
axis([-200 200 y_lim(1) y_lim(2)])
str = ['\it{Mean angle: }' '\bf ' num2str(ang_ibsang)];
text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','red')
str = ['\it{Mean resultant length: }' '\bf ' num2str(mvl_ibsang)];
text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/12,str,'Color','red')
str = ['\it{n: }' '\bf ' num2str(n_ibsang)];
text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/6,str,'Color','red')
subplot(2,2,3)      % cycle first single spikes
bar(cnts,nm_sspoang'/length(aang_sspoang))
title(gca,[titlestr ' cycle first single spikes'])
y_lim = ylim;
axis([-200 200 y_lim(1) y_lim(2)])
str = ['\it{Mean angle: }' '\bf ' num2str(ang_sspoang)];
text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','red')
str = ['\it{Mean resultant length: }' '\bf ' num2str(mvl_sspoang)];
text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/12,str,'Color','red')
str = ['\it{n: }' '\bf ' num2str(n_sspoang)];
text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/6,str,'Color','red')
subplot(2,2,4)
bar(xout_cycnb,nm_cycnb)
title(gca,[titlestr ' distr. of burst no./cycle'])
y_lim = ylim;
str = ['\it{Cycles missed: }' '\bf ' num2str(length(find(cycnb==0))/length(cycnb))];
text(2,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/6,str,'Color','red')

figure(H4)      % STA
set(gcf,'Position',[1 31 1278 920])
subplot(2,2,1)      % all spikes
time = linspace(-wn/sr/2,wn/sr/2,length(sta_as));
plot(time,sta_as,'LineWidth',1.5)
title(gca,[titlestr ' all spikes'])
y_lim = ylim;
str = ['\it{Max-mean: }' '\bf ' num2str(sta_index1_as)];
text(0.1,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','red')
str = ['\it{Max: }' '\bf ' num2str(sta_index2_as)];
text(0.1,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/12,str,'Color','red')
str = ['\it{n: }' '\bf ' num2str(nn_as)];
text(0.1,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/6,str,'Color','red')
subplot(2,2,2)      % all first spikes
time = linspace(-wn/sr/2,wn/sr/2,length(sta_afsp));
plot(time,sta_afsp,'LineWidth',1.5)
title(gca,[titlestr ' all first spikes'])
y_lim = ylim;
str = ['\it{Max-mean: }' '\bf ' num2str(sta_index1_afsp)];
text(0.1,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','red')
str = ['\it{Max: }' '\bf ' num2str(sta_index2_afsp)];
text(0.1,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/12,str,'Color','red')
str = ['\it{n: }' '\bf ' num2str(nn_afsp)];
text(0.1,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/6,str,'Color','red')
subplot(2,2,3)      % burst first spikes
time = linspace(-wn/sr/2,wn/sr/2,length(sta_fs));
plot(time,sta_fs,'LineWidth',1.5)
title(gca,[titlestr ' burst first spikes'])
y_lim = ylim;
str = ['\it{Max-mean: }' '\bf ' num2str(sta_index1_fs)];
text(0.1,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','red')
str = ['\it{Max: }' '\bf ' num2str(sta_index2_fs)];
text(0.1,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/12,str,'Color','red')
str = ['\it{n: }' '\bf ' num2str(nn_fs)];
text(0.1,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/6,str,'Color','red')
subplot(2,2,4)      % single spikes
time = linspace(-wn/sr/2,wn/sr/2,length(sta_sp));
plot(time,sta_sp,'LineWidth',1.5)
title(gca,[titlestr ' single spikes'])
y_lim = ylim;
str = ['\it{Max-mean: }' '\bf ' num2str(sta_index1_sp)];
text(0.1,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','red')
str = ['\it{Max: }' '\bf ' num2str(sta_index2_sp)];
text(0.1,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/12,str,'Color','red')
str = ['\it{n: }' '\bf ' num2str(nn_sp)];
text(0.1,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/6,str,'Color','red')

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
fns = [fname(1:end-4) '_BURSTVSPHASE1.fig'];
saveas(H1,fns)
fns = [fname(1:end-4) '_BURSTVSPHASE1.jpg'];
saveas(H1,fns)
fns = [fname(1:end-4) '_BURSTVSPHASE2.fig'];
saveas(H2,fns)
fns = [fname(1:end-4) '_BURSTVSPHASE2.jpg'];
saveas(H2,fns)
fns = [fname(1:end-4) '_PHASEHIST.fig'];
saveas(H3,fns)
fns = [fname(1:end-4) '_PHASEHIST.jpg'];
saveas(H3,fns)
fns = [fname(1:end-4) '_STA.fig'];
saveas(H4,fns)
fns = [fname(1:end-4) '_STA.jpg'];
saveas(H4,fns)
fns = [fname(1:end-4) '_VMFIT1.fig'];
saveas(H5,fns)
fns = [fname(1:end-4) '_VMFIT1.jpg'];
saveas(H5,fns)
fns = [fname(1:end-4) '_VMFIT2.fig'];
saveas(H6,fns)
fns = [fname(1:end-4) '_VMFIT2.jpg'];
saveas(H6,fns)
fns = [fname(1:end-4) '_CYCLEFIRST.fig'];
saveas(H8,fns)
fns = [fname(1:end-4) '_CYCLEFIRST.jpg'];
saveas(H8,fns)
fn = [fname(1:end-4) '_STA.mat'];
save(fn,'sta_fs','sta_sp','sta_as','sta_afsp','sta_index1_fs','sta_index1_sp',...
    'sta_index1_as','sta_index1_afsp','sta_index2_fs','sta_index2_sp',...
    'sta_index2_as','sta_index2_afsp');
fn = [fname(1:end-4) '_PHASE.mat'];
save(fn,'aang_fs','aang_sp','aang_as','aang_afsp','aang_ibsang',...
    'aang_sspoang','aang_allang','cycnb','r_fs','r_sp','r_as','r_afsp');
fn = [fname(1:end-4) '_BURSTPHASE.mat'];
save(fn,'H1ibspno','H1ibfr','H1bl','H1aang_fs');
close all

cd(resdir2)
fn = [fname(1:end-4) '_CLUST2.mat'];
save(fn,'Burstiness','IntraBurstFrequency','IntraBurstSpikeNumber','BurstLength',...
    'BurstFrequency','IsiMatrix','FiringRate','Burst','vdisc')

cmps = strread(fname(1:end-4),'%s','delimiter','_');
xlsname = cmps{1};
if exist([xlsname '.xls'],'file')
    ntx = xlsread(xlsname,'sheet1');
    nty = xlsread(xlsname,'sheet2');
    [ntz mtz atz] = xlsread(xlsname,'sheet3');
else
    ntx = [];
    nty = [];
    ntz = [];
    atz = {};
    hrow = {[] 'Bness' 'IBF mean' 'IBF sd' 'IBSN mean' 'IBSN sd' 'BL mean' 'BL sd' 'BF'};
    hrowy = {[] 'IBF' 'IBSN' 'BL'};
    xlswrite(xlsname,hrow,'sheet1','A1')
    xlswrite(xlsname,hrowy,'sheet2','A1')
end
mt = {fname};
ntx2 = [Burstiness IntraBurstFrequency.mean IntraBurstFrequency.sd...
    IntraBurstSpikeNumber.mean IntraBurstSpikeNumber.sd...
    BurstLength.mean BurstLength.sd BurstFrequency];
nty2 = [IntraBurstFrequency.all' IntraBurstSpikeNumber.all' BurstLength.all'];
ntz2 = [IsiMatrix.mean IsiMatrix.sd IsiMatrix.num];
hrowz = cell(1,1+2*size(IsiMatrix.mean,2)+size(IsiMatrix.num,2));
hrowz{2} = 'mean';
hrowz{1+size(IsiMatrix.mean,2)+1} = 'sd';
hrowz{1+2*size(IsiMatrix.mean,2)+1} = 'num';
str = ['A' num2str(size(ntx,1)+2)];
xlswrite(xlsname,mt,'sheet1',str)
str = ['B' num2str(size(ntx,1)+2)];
xlswrite(xlsname,ntx2,'sheet1',str)
str = ['A' num2str(size(nty,1)+2)];
xlswrite(xlsname,mt,'sheet2',str)
str = ['B' num2str(size(nty,1)+2)];
xlswrite(xlsname,nty2,'sheet2',str)
str = ['A' num2str(size(atz,1)+1)];
xlswrite(xlsname,hrowz,'sheet3',str)
str = ['A' num2str(size(atz,1)+2)];
xlswrite(xlsname,mt,'sheet3',str)
str = ['B' num2str(size(atz,1)+2)];
xlswrite(xlsname,ntz2,'sheet3',str)



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

% -------------------------------------------------------------------------
function [ang_ibsang ang_sspoang ang_allang cycnb] = laphaseb(feeg,ahee,sspo,ibs,vb1,sr)
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