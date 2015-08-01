function aall
%AALL   Time-resolved phase histogram, STA, burst statistics and entropy.
%   AALL calculates phase histograms, Spike Triggered Average, burst
%   parameters and entropy measures for half-second-long time windows.
%   Results are saved in the directories specified in the first part of the
%   program code.
%
%   Output:
%       1. Burstiness, intraburst frequency, intraburst spike number, burst
%       frequency and firing rate are displayed and saved in a multiple
%       plot. A mat and an excel file containing the variables are also 
%       saved by the program.
%       2. Color-coded phase histograms, mean angle and mean resultant
%       length are displayed and saved in a multiple plot.
%       3. Phase-specific firing rate for 60 degree phase segments are
%       plotted and saved.
%       4. Color-coded Spike Triggered Average amplitude (with 2-sec. 
%       window size) and the two types of STA indices (see ASTANORM for
%       details) are plotted and saved.
%       5. Unit and EEG wavelet entropy, mutual information, EEG->unit and
%       unit->EEG uncertainty coefficients are saved in mat and fig files.
%       Comparisons between real data and shuffled segment controls (see
%       ENTRYRUN3C_CONT2 for details) with Wilcoxon signedrank test are
%       saved in boxplot figures.
%
%   See also ENTRYRUN3C_CONT2, AENTRY, APHASE2, ASTANORM and ACLUSTRUN.

% Directories
global DATAPATH
inpdir1 = 'X:\In_Vivo\_analysis\acsady_csoport\auj_istvan\mat_ket_xyl_3\';
inpdir2 = [DATAPATH 'Andi\Disc\ketxyl\'];
inpdir3 = [DATAPATH 'Andi\Ketxyl\Cluster\mat2\'];
resdir1 = [DATAPATH 'Andi\Ketxyl\BurSTat\mat\'];
resdir2 = [DATAPATH 'Andi\Ketxyl\BurSTat\fig\'];
resdir3 = [DATAPATH 'Andi\Ketxyl\Length\'];
resdir4 = [DATAPATH 'Andi\Ketxyl\BurSTat\xls\'];
resdir5 = [DATAPATH 'Andi\Ketxyl\PhaseSTat\Phase\'];
resdir6 = [DATAPATH 'Andi\Ketxyl\PhaseSTat\STA\'];
resdir7 = [DATAPATH 'Andi\Ketxyl\EntrySTat_gamma\'];
mm = pwd;
cd(inpdir1)

% Filelist
[files1 files_short1] = filelist(inpdir1);
[files2 files_short2] = filelist2(inpdir2);
files_short = intersect(files_short1,files_short2);
sf = length(files_short);

% Progress indicator
[wb,awb1,awb2] = waitbar2([0 0],'Running AALL...');
global WB
WB(end+1) = wb;

% MAIN
sr = 20000;         % sapling rate
dsr = 1000;         % downsampling for phase hist., STA and entropy calculation
nextcol = 1;
edges = -180:20:180;     % edges for phase histogram
wn = 2 * dsr;    % 2 sec. window for STA
for o = 1:sf
    fname = files_short{o};
    cmps = strread(fname,'%s','delimiter','_');
    if length(cmps) < 3
        strw = [cmps{1} ' ' cmps{2}];
    else
        strw = [cmps{1} ' ' cmps{2} ' ' cmps{3}];
    end
    waitbar2([(o-1)/sf 0],wb,strw);
    ff3 = [inpdir3 fname(1:end-4) '_CLUST2'];
    try
        load(ff3)       % load burst analysis results
    catch
        continue
    end
    ff2 = [inpdir2 fname(1:end-4) '_d.mat'];
    load(ff2)           % load discriminated unit
    if length(vdisc) < 50
        continue
    end
    ff = [inpdir1 fname];
    load(ff)            % load original data
    len = length(data);
    eeg = data(:,2)';
    clear data
    
%     efflen = (vdisc(end) - vdisc(1)) / sr;
%     frate = length(vdisc) / efflen;
%     seglen = 50 / frate;        % adaptive segment length: should contain at least 50 spikes on average
    seglen = 30 * sr;        % 30 sec. long segments
    lenr = floor(len/seglen);       % preallocation
    ind1 = [1:seglen:len];
    ind2 = ind1 + seglen -1;
    burstiness = zeros(1,lenr);
    burstfreq = zeros(1,lenr);
    frate = zeros(1,lenr);
    vburst = vdisc(Burst);
    ang = zeros(1,lenr);
    mvl = zeros(1,lenr);
    Phist = zeros(length(edges)-1,lenr);
    phfr = zeros(6,lenr);
    sta_index1 = zeros(1,lenr);
    sta_index2 = zeros(1,lenr);
    Sta = zeros(wn+1,lenr);
    rHx = [];
    rHy = [];
    rHxy = [];
    rIxy = [];
    rUxy = [];
    rUyx = [];
    cHx = [];
    cHy = [];
    cHxy = [];
    cIxy = [];
    cUxy = [];
    cUyx = [];
    for k = 1:lenr
        vd = vdisc(vdisc>ind1(k)&vdisc<ind2(k)) - ind1(k);
        loceeg = eeg(ind1(k):ind2(k));
        
% Burst statistics
        lvb = locvburst(vburst,ind1(k),ind2(k));
        burstnum = size(lvb,2);
        intraburstiv = [];
        intraburstnum = zeros(1,burstnum);
        for j = 1:burstnum      % computing intraburstiv
            b = vdisc(vdisc>=lvb(1,j)&vdisc<=lvb(2,j));
            db = diff(b);
            intraburstiv = [intraburstiv db];
            intraburstnum(j) = length(b);   % intraburst spike number
        end
        burstiness(k) = (length(intraburstiv) + burstnum) / length(vd);
        burstlength = (lvb(2,:) - lvb(1,:)) / sr;
        if ~isempty(intraburstnum)
            intraburstfreq(k) = mean((intraburstnum-1)./burstlength);
            ibspno(k) = mean(intraburstnum);
        else
            intraburstfreq(k) = NaN;
            ibspno(k) = NaN;
        end
        burstfreq(k) = 20000 * (burstnum - 1)  / (vdisc(Burst(2,end)) - vdisc(Burst(1,1)));
        efflen = (vd(end) - vd(1)) / sr;
        frate(k) = length(vd) / efflen;
        
% Phase histograms
        eeg2 = loceeg(1:20:end);    % downsample on 1000 Hz
        vdisc2 = round(vd/20);
        dsr = 1000;

        [aang cyclen1 cyclen2 cl phlen] = aphase3(eeg2,vdisc2,dsr);    % PHASE
        n = length(aang);
        ftm = sum(exp(1).^(i*aang)) / n;    % first trigonometric moment
        ang(k) = angle(ftm) * 180 / pi;   % mean angle
        mvl(k) = abs(ftm);     % mean resultant length
        aang = aang * 180 / pi;
        
        nm = histc(aang,edges);   % phase histogram
        Phist(:,k) = nm(1:end-1)';
        
        grnm = [sum(nm(1:3)) sum(nm(4:6)) sum(nm(7:9)) sum(nm(10:12)) sum(nm(13:15)) sum(nm(16:18))];
        phfr(1:6,k) = grnm ./ phlen;

% Spike Triggered Average
        [sta sta_index1(k) sta_index2(k)] = astanorm(vdisc2,eeg2,wn);    % STA
        Sta(:,k) = sta';
        
% Entropy
%         [rHxt,rHyt,rHxyt,rIxyt,rUxyt,rUyxt,cHxt,cHyt,cHxyt,cIxyt,cUxyt,cUyxt] = aentry2(eeg2,vd,sr);
%         rHx = [rHx rHxt];
%         rHy = [rHy rHyt];
%         rHxy = [rHxy rHxyt];
%         rIxy = [rIxy rIxyt];
%         rUxy = [rUxy rUxyt];
%         rUyx = [rUyx rUyxt];
%         cHx = [cHx cHxt];
%         cHy = [cHy cHyt];
%         cHxy = [cHxy cHxyt];
%         cIxy = [cIxy cIxyt];
%         cUxy = [cUxy cUxyt];
%         cUyx = [cUyx cUyxt];
        
        waitbar2([(o-1)/sf k/lenr],wb,strw);
    end
    
% Save
    figure;     % Burst statistics - burstiness
    subplot(5,1,1)
    plot(burstiness);
    set(gca,'YLim',[0 1])
    title('Burstiness')
    
    subplot(5,1,2)      % intraburst frequency
    plot(intraburstfreq);
    set(gca,'YLim',[100 600])
    title('Intraburst Frequency')
    
    subplot(5,1,3)      % intraburst spike number
    plot(ibspno);
    set(gca,'YLim',[0 10])
    title('Intraburst Spike Number')
    
    subplot(5,1,4)      % burst fequency
    plot(burstfreq);
    set(gca,'YLim',[0 0.1])
    title('Burst Frequency')
    
    subplot(5,1,5)      % firing rate
    plot(frate);
    set(gca,'YLim',[0 10])
    title('Firing Rate')
    
    fn = [resdir2 fname(1:end-4) '_BURST'];     % save figure
    saveas(gcf,fn)
    
    fn = [resdir1 fname(1:end-4) '_BURSTAT'];   % save mat file
    save(fn,'burstiness','intraburstfreq','ibspno','burstfreq','frate','seglen')
    
    fn = [resdir3 fname(1:end-4) '_LEN'];       % save data length
    save(fn,'len')
    
    colref1 = [colno2colref(nextcol) '1'];      % save excel file
    colref2 = [colno2colref(nextcol) '2'];
    xlswrite([resdir4 'burstiness.xls'],{fname},'sheet1',colref1)
    xlswrite([resdir4 'burstiness.xls'],burstiness','sheet1',colref2)
    xlswrite([resdir4 'intraburstfreq.xls'],{fname},'sheet1',colref1)
    xlswrite([resdir4 'intraburstfreq.xls'],intraburstfreq','sheet1',colref2)
    xlswrite([resdir4 'ibspno.xls'],{fname},'sheet1',colref1)
    xlswrite([resdir4 'ibspno.xls'],ibspno','sheet1',colref2)
    xlswrite([resdir4 'burstfreq.xls'],{fname},'sheet1',colref1)
    xlswrite([resdir4 'burstfreq.xls'],burstfreq','sheet1',colref2)
    xlswrite([resdir4 'frate.xls'],{fname},'sheet1',colref1)
    xlswrite([resdir4 'frate.xls'],frate','sheet1',colref2)
    nextcol = nextcol + 1;
    
    figure;     % Phase histograms
    subplot(3,1,1)
    imagesc(Phist)
    subplot(3,1,2)
    plot(ang)
    subplot(3,1,3)
    plot(mvl)
    set(gca,'YLim',[0 1])
    fn = [resdir5 fname(1:end-4) '_PHASE'];
    saveas(gcf,fn)
    
    figure;     % Phase firing rate
    for t = 1:6
        subplot(6,1,t)
        plot(phfr(t,:))
        set(gca,'YLim',[0 15])
        llim = -180 + (t - 1) * 60;
        ulim = llim + 60;
        tstr = [num2str(llim) ' - ' num2str(ulim)];
        title(tstr)
    end
    fn = [resdir5 fname(1:end-4) '_PHFR'];
    saveas(gcf,fn)
    
    
    figure;     % STA
    subplot(3,1,1)
    imagesc(Sta)
    subplot(3,1,2)
    plot(sta_index1)
    subplot(3,1,3)
    plot(sta_index2)
    fn = [resdir6 fname(1:end-4) '_STA'];
    saveas(gcf,fn)
%     
%     fn = [resdir7 fname(1:end-4) '_ENTRY'];         % Entropy - save mat file
%     save(fn,'rHx','rHy','rHxy','rIxy','rUxy','rUyx','cHx','cHy','cHxy','cIxy','cUxy','cUyx');
% 
%     x1 = 0;     % uncert. coeff.
%     x2 = len / sr;
%     xa = linspace(x1,x2,length(rHx));
%     figure
%     subplot(2,1,1)
%     L(1) = plot(xa,rUxy);
%     set(gca,'YLim',[0 0.5])
%     hold on
%     L(2) = plot(xa,rUyx,'r');
%     set(gca,'YLim',[0 0.5])
%     y_lim = ylim;
%     axis([x1 x2 y_lim(1) y_lim(2)]);
%     legend(L,'eeg->unit','unit->eeg');
%     subplot(2,1,2)
%     plot(xa,cUxy)
%     set(gca,'YLim',[0 0.5])
%     hold on
%     plot(xa,cUyx,'r')
%     set(gca,'YLim',[0 0.5])
%     y_lim = ylim;
%     axis([x1 x2 y_lim(1) y_lim(2)]);
%     
%     fn = [resdir7 fname(1:end-4) '_UNCERT'];
%     saveas(gcf,fn)
%     
%     figure;     % mutual information
%     subplot(2,1,1)
%     plot(xa,rIxy)
%     set(gca,'YLim',[0 2])
%     y_lim = ylim;
%     axis([x1 x2 y_lim(1) y_lim(2)]);
%     subplot(2,1,2)
%     plot(xa,cIxy)
%     set(gca,'YLim',[0 2])
%     y_lim = ylim;
%     axis([x1 x2 y_lim(1) y_lim(2)]);
%     
%     fn = [resdir7 fname(1:end-4) '_MUTINF'];
%     saveas(gcf,fn)
%     
%     figure;     % mutual information boxplot and Wilcoxon signedrank test
%     boxplot([rIxy cIxy],[zeros(size(rIxy)) ones(size(cIxy))],'labels',...
%         [{'mutual information'} {'control'}])
%     [Wp_eu,Wh_eu] = b_signrank2(rIxy,cIxy,'alpha',0.05);
%     if Wh_eu
%         clr = 'red';
%     else
%         clr = 'black';
%     end
%     y_lim = ylim;
%     x_lim = xlim;
%     tpos1 = x_lim(1) + (x_lim(2) - x_lim(1)) / 2;
%     tpos2 = y_lim(1) + (y_lim(2)-y_lim(1)) * 4 / 5;
%     text(tpos1,tpos2,num2str(Wp_eu),'Color',clr,'Horizontalalignment','center')
%     
%     fn = [resdir7 fname(1:end-4) '_BOXmutinf'];
%     saveas(gcf,fn)
%         
%     figure;     % entropy boxplot
%     boxplot([rHx cHx rHy cHy],[zeros(size(rHx)) ones(size(cHx))...
%         2*ones(size(rHy)) 3*ones(size(cHy))],'labels',[{'unit entropy'}...
%         {'control'} {'EEG entropy'} {'control'}])
%     
%     fn = [resdir7 fname(1:end-4) '_BOXentropy'];
%     saveas(gcf,fn)
%     
%     figure;     % uncertainty coeff. boxplot and Wilcoxon signedrank test
%     boxplot([rUxy cUxy rUyx cUyx],[zeros(size(rUxy)) ones(size(cUxy))...
%         2*ones(size(rUyx)) 3*ones(size(cUyx))],'labels',[{'EEG->unit'}...
%         {'control'} {'unit->EEG'} {'control'}])
%     [Wp_eu,Wh_eu] = b_signrank2(rUxy,cUxy,'alpha',0.05);
%     if Wh_eu
%         clr = 'red';
%     else
%         clr = 'black';
%     end
%     y_lim = ylim;
%     x_lim = xlim;
%     tpos1 = x_lim(1) + (x_lim(2) - x_lim(1)) / 4;
%     tpos2 = y_lim(1) + (y_lim(2)-y_lim(1)) * 4 / 5;
%     text(tpos1,tpos2,num2str(Wp_eu),'Color',clr,'Horizontalalignment','center')
%     [Wp_eu,Wh_eu] = b_signrank2(rUyx,cUyx,'alpha',0.05);
%     if Wh_eu
%         clr = 'red';
%     else
%         clr = 'black';
%     end
%     y_lim = ylim;
%     x_lim = xlim;
%     tpos1 = x_lim(1) + (x_lim(2) - x_lim(1)) * 3 / 4;
%     tpos2 = y_lim(1) + (y_lim(2)-y_lim(1)) * 4 / 5;
%     text(tpos1,tpos2,num2str(Wp_eu),'Color',clr,'Horizontalalignment','center')
%     
%     fn = [resdir7 fname(1:end-4) '_BOXuncert'];
%     saveas(gcf,fn)
    
    close all
    
end
cd(mm)
close(wb)



% -------------------------------------------------------------------------
function [files2 files2_short] = filelist(inpdir)

% List of filenames
files = dir(inpdir);
files = files(3:end);
files2 = struct('name',[],'date',[],'bytes',[],'isdir',[],'datenum',[]);
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
files2 = struct('name',[],'date',[],'bytes',[],'isdir',[],'datenum',[]);
files2_short = {};
for i = 1:length(files)
    if ~files(i).isdir
        files2(end+1) = files(i);
        files2_short{end+1} = [files(i).name(1:end-6) '.mat'];
    end
end
files2 = files2(2:end);



% -------------------------------------------------------------------------
function vburst = locvburst(vburst,ind1,ind2)

VBurst = vburst;
vburst = VBurst .* (VBurst>ind1) .* (VBurst<ind2);
if isequal(size(VBurst),[1,2])
    VBurst = VBurst';
end
vburst = VBurst .* (VBurst>ind1) .* (VBurst<ind2);  % restrict vburst to the window
if isempty(vburst)
    return
end
vb1 = [vburst(1,:)>0 vburst(2,:)>0] .* [[1:size(vburst,2)] (-1)*[1:size(vburst,2)]];
if sum(vb1) > 0     % handle the case when burst exceeds the window
    ind = find(vburst,1,'last')+1;
    vburst(ind) = ind2;
    vb1 = [vburst(1,:)>0 vburst(2,:)>0] .* [[1:size(vburst,2)] (-1)*[1:size(vburst,2)]];
end
if sum(vb1) < 0
    ind = find(vburst,1,'first')-1;
    vburst(ind) = ind1;
end
fV1 = find(VBurst>ind2,1,'first');
fV2 = find(VBurst<ind1,1,'last');
if fV1 - fV2 == 1 & mod(fV1,2) == 0
    vburst = [ind1 ind2];   % when one burst exceeds the window on both sides
end
vburst = vburst(find(vburst));
vburst = reshape(vburst,2,length(vburst)/2);



% -------------------------------------------------------------------------
function [ang cyclen1 cyclen2 cl] = aphase2(eeg,vdisc,sr)
%APHASE2    Phase angles for unit relative to EEG.
%   A,C1,C2,CL = APHASE2(EEG,VDISC,SR) calculates Hilbert phase angles (A)
%   for discriminated unit (VDISC) relative to EEG, when sampling frequency
%   is given in SR. Cicles not fulfilling the following 2 criteria are
%   discarded: (i) EEG amp. higher then 2SD; (ii) min. 100 ms length (half
%   wavelength of filter cutoff freq.). Mean cycle length before and after
%   applying the criteria are returned in C1 and C2. CL cell contains length
%   of the following cycles:
%       1. original
%       2. discarded upon 1st crit.
%       3. discarded upon 2nd crit.
%       4. all discarded
%       5. all remaining.
%
%   See also HILBERT.

% Filtering EEG
nqf = sr / 2;
flt = fir1(4096,5/nqf,'low');      % lowpass filtering on 5 Hz
feeg = filtfilt(flt,1,eeg);

% Hilbert transformation
ahee = angle(hilbert(feeg));

% Check SWS criteria:
% 1. discard cicles with EEG amp. lower then 2SD
% 2. discard cicles shorter then 100 ms (half wavelength of filter cutoff freq.)
fn = find(-diff(ahee)>2*pi-0.1);
cyclen1 = mean(diff(fn)) / sr * 1000;   % cycle length in ms
sd = std(feeg);
inx = find(vdisc<fn(1));
cl1 = [];
cl4 = [];
cl5 = [];
cl6 = [];
for k = 1:length(fn)-1
    seeg = feeg(fn(k):fn(k+1));
    axs = max(seeg) - min(seeg);
    sahee = ahee(fn(k):fn(k+1));
    if (axs < 2 * sd)  | (fn(k+1) - fn(k) < 0.25 * sr)
        inx = [inx find(vdisc>fn(k)&vdisc<fn(k+1))];
        cl5(end+1) = (fn(k+1) - fn(k)) / sr * 1000;   % discarded cycles' length in ms;
    else
        cl6(end+1) = (fn(k+1) - fn(k)) / sr * 1000;   % remaining cycles' length in ms;
    end
    if axs < 2 * sd
        cl1(end+1) = (fn(k+1) - fn(k)) / sr * 1000;   % cycle length in ms;
    end
    if fn(k+1) - fn(k) < 0.25 * sr
        cl4(end+1) = (fn(k+1) - fn(k)) / sr * 1000;   % cycle length in ms;
    end
end
inx = [inx find(vdisc>fn(end))];
vdisc(inx) = [];
cyclen2 = mean(cl6) / sr * 1000;   % cycle length in ms
cl{1} = diff(fn) ./ sr .* 1000;   % original cycle length in ms
cl{2} = cl1;    % discarded upon 1st crit.
cl{3} = cl4;    % discarded upon 2nd crit.
cl{4} = cl5;    % all discarded
cl{5} = cl6;    % all remaining
ang = ahee(vdisc);