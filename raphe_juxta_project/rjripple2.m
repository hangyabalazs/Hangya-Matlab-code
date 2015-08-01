function rjripple2
%RJRIPPLE2   Firing of MRN cells during ripples.
%   RJRIPPLE2 calulates firing histograms, raw unit overlays and filtered
%   (90-140 Hz) LFP overlays centered on ripple maxima. It calls SSPW for
%   ripple detection. It also calculates normalized ripple histograms (see
%   Klausberger et al., 2003).
%
%   See also SSPW and RJRIPPLE.

% Directories
global DATADIR
global DATAPATH
inpdir = [DATADIR 'raphe_matfiles\raphe_juxta_files_discriminated\'];
inpdir2 = [DATADIR 'raphe_matfiles\raphe_juxta_files\'];
thetadir = [DATAPATH 'Raphe\raphe_juxta\Wavelet\theta_segments\'];
nothdir = [DATAPATH 'Raphe\raphe_juxta\Wavelet\nontheta_segments\'];
sharpdir = [DATAPATH 'Raphe\raphe_juxta\Wavelet\sharpwave_segments\'];
resdir = [DATAPATH 'Raphe\raphe_juxta\Ripple2\'];
mm = pwd;
cd(resdir)

% Filelist
[files files_short] = b_filelist(inpdir);
sf = length(files_short);

% Progress indicator
wb = waitbar(0,'Running RJRIPPLE...','Position',[360 250 275 50]);
global WB
WB(end+1) = wb;

% Main
sr = 10000;     % sampling rate
nqf = sr / 2;   % Nyquist frequency
for o = 12:sf
    fname = files(o).name;
    cmps = strread(fname,'%s','delimiter','_');
    titlestr = [cmps{1} ' ' cmps{2}];
    ff = [inpdir fname];
    load(ff)
    ff2 = [inpdir2 fname(1:end-6) '.mat'];
    load(ff2)
    unit = Unit.values;
    eeg_times = linspace(1,length(eeg)/sr,length(eeg));
    
    [bg fin feeg] = sspw(eeg,sr);
    if isempty(bg)
        continue
    end
    spwno = length(bg);
    pks = zeros(1,spwno);
    figure
    wn = 10000;
    z = zeros(1,2*wn);
    for k = 1:spwno
        lbig = bg(k);
        lfin = fin(k);
        if lfin > length(feeg)
            continue
        end
        lfeeg = feeg(lbig:lfin);
        pk = find(lfeeg==max(lfeeg));
        cnt = lbig + pk(1);
        pks(k) = cnt;
        if cnt - wn < 1 || cnt + wn > length(feeg)
            continue
        end
        loceeg = feeg(cnt-wn:cnt+wn);
        S1 = subplot(3,1,1);
        hold on
        plot(loceeg)
        S2 = subplot(3,1,2);
        hold on
        locunit = unit(cnt-wn:cnt+wn);
        plot(locunit)
        locvdisc = (vdisc(vdisc>cnt-wn&vdisc<cnt+wn)) - cnt + wn;
        z(locvdisc) = 1;
    end
    linkaxes([S1 S2],'x')
    z2 = reshape(z,2*wn/20,20);
    subplot(3,1,3)
    bar(1:20,sum(z2))
    title(titlestr)
    fn = [fname(1:end-4) '_RIPPLE'];
    saveas(gcf,fn);
    
% Normalized ripple histogram
    term = spwno;   % terminate if spiking ceased
    pterm = 1;      % begin when spiking starts
    for k = 1:spwno
        be = bg(k);
        vdpre = vdisc(find(vdisc<be,1,'first'));
        if isempty(vdpre)
            pterm = k + 1;
        end
        fe = fin(k);
        vdpost = vdisc(find(vdisc>fe,1,'first'));
        if isempty(vdpost)
            term = k - 1;
            break
        end
    end
    nmn = zeros(term-pterm+1,32);
    nrp = 1000;
    slev = 0.05 * nrp;  % 5% significance level
    prnmn = zeros(term-pterm+1,nrp,32);
    for k = pterm:term
        iv1 = (pks(k) - bg(k)) / 4;
        iv2 = (fin(k) - pks(k)) / 4;
        ivm = (iv1 + iv2) / 2;
        w1 = bg(k) - 12 * ivm;
        w2 = fin(k) + 12 * ivm;
        edges = [w1:ivm:(bg(k)-ivm) bg(k):iv1:(pks(k)-iv1)...
            pks(k):iv2:(fin(k)-iv2) fin(k):ivm:w2];
        [nm,xout] = histc(vdisc,edges);
        nm = nm(1:end-1);
        nmn(k,:) = nm ./ [ones(1,12)*ivm ones(1,4)*iv1 ones(1,4)*iv2 ones(1,12)*ivm] * sr; % in spikes/sec
        
        vdpre = vdisc(find(vdisc<edges(1),1,'last'));    % permutation test
        vdpost = vdisc(find(vdisc>edges(end),1,'first'));
        lvd = vdisc(vdisc>=vdpre&vdisc<=vdpost);
        dlvd = diff(lvd);
        lisi = length(dlvd);
        for kk = 1:nrp
            rp = randperm(lisi);
            vdt = vdpre + cumsum(dlvd(rp));
            rnm = histc(vdt,edges)';
            rnm = rnm(1:end-1);
            rnm = rnm(:)';
            prnmn(k,kk,:) = rnm ./ [ones(1,12)*ivm ones(1,4)*iv1 ones(1,4)*iv2 ones(1,12)*iv2] * sr;
        end
        
    end
    mprnmn = squeeze(mean(prnmn,1));
    mnmn = mean(nmn);
    prpeak = mprnmn(:,13:20);
    rpeak = mean(prpeak,2);
    plpeak = mnmn(:,13:20);
    lpeak = mean(plpeak,2);
    nx = find(sort(rpeak,'ascend')<lpeak,1,'last');
    str = ' ';
    if nx > 1000 - slev
        str = [str ' single peak'];
    elseif nx < slev
        str = [str ' anti-sharp-wave'];
    end
    prpeak1 = mprnmn(:,13:16);
    prpeak2 = mprnmn(:,17:20);
    rpeak1 = mean(prpeak1,2);
    rpeak2 = mean(prpeak2,2);
    rdiffp = rpeak1 - rpeak2;
    plpeak1 = mnmn(:,13:16);
    plpeak2 = mnmn(:,17:20);
    ldiffp = mean(plpeak1) - mean(plpeak2);
    nx = find(sort(rdiffp,'ascend')<ldiffp,1,'last');
    if nx > 1000 - slev
        str = [str ' biphasic'];
    end
    
    figure
    bar(1:32,mnmn,'FaceColor','black');
    hold on
    bar(13:16,mnmn(13:16),'FaceColor','red');
    bar(17:20,mnmn(17:20),'FaceColor','blue');
    bar(21:32,mnmn(21:32),'FaceColor','black');
    x_lim = xlim;
    y_lim = ylim;
    text(x_lim(1)+(x_lim(2)-x_lim(1))*0.7,y_lim(1)+(y_lim(2)-y_lim(1))*0.7,str)
    ts = fname(1:end-4);
    ts(ts=='_') = ' ';
    title(ts)
    fns = [fname(1:end-4) '_RIPPLEk.fig'];  % save (k after Klausberger)
    saveas(gcf,fns)
    
    close all
    waitbar(o/sf)
end
close(wb)
cd(mm)