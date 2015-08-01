function czripple
%CZRIPPLE   Firing of interneurons during ripples.
%   CZRIPPLE calculates firing histograms during ripples centered on (i)
%   ripple maximum, (ii) ripple start, (iii) ripple end. Normalized ripple
%   histograms (see Klausberger et al., 2003) are also calculated.
%   Significant activation, inactivation or biphasic response is tested
%   using a permutation test.
%
%   Reference:
%   Klausberger T, Magill PJ, Márton LF, Roberts JD, Cobden PM, Buzsáki G,
%   Somogyi P (2003) Brain-state- and cell-type-specific firing of
%   hippocampal interneurons in vivo. Nature  421:844-848.
%
%   See also CZPHASE, RAPHEPREVIEW and SSPW.

% Directories
global DATAPATH
inpdir_eeg1 = [DATAPATH 'Czurko\czurko_EEG\'];
inpdir_eeg2 = [DATAPATH 'Czurko\czurko_EEG\resample\'];
inpdir_unit = [DATAPATH 'Czurko\discriminated2\'];
resdir = [DATAPATH 'Czurko\czurko_EEG\neg4\'];
mm = pwd;
cd(resdir)

% Load
xlsname = [inpdir_eeg1 'EEG2.xls'];
[ntz mtz] = xlsread(xlsname,'neg');
sf = size(mtz,1);   % number of pairs
for o = 1:1
    dr = [inpdir_unit mtz{o,4} '\' mtz{o,1} '\'];
    fn = [dr mtz{o,2} '.mat'];
    load(fn)        % load unit
    vdisc = data;
    
    fn = [inpdir_eeg2 'EEG_' mtz{o,3} '_' mtz{o,1} '_rs.mat'];
    load(fn)        % load resampled EEG (sampled on 10 kHz)
    fn = [inpdir_eeg1 'EEG_' mtz{o,3} '_' mtz{o,1} '.mat'];
    load(fn)        % load original EEG
    eval(['Eeg = ' mtz{o,3} ';']);
    eval(['clear ' mtz{o,3}]);
%     eeg = Eeg.values;
%     sr = 1 / Eeg.interval;
    sr = 10000;
    eeg_start = Eeg.start;
    eeg_end = eeg_start + (length(eeg) - 1) / sr;
    eeg_times = eeg_start:1/sr:eeg_end;
    
% Downsample EEG
    eeg2 = eeg(1:5:end);
    eeg_times2 = eeg_times(1:5:end);
    sr2 = sr / 5;
        
% Ripples
    [bg fin feeg] = sspw(eeg,sr);
    figure
    plot(eeg)
    spwno = length(bg);
    for k = 1:spwno
        hold on
        plot([bg(k):fin(k)],eeg(bg(k):fin(k)),'r')
    end
    ts = mtz{o,1};
    ts(ts=='_') = ' ';
    fns = [ts '_RIPPLES.fig'];  % save
    saveas(gcf,fns)
    
% Ripples - centered on peak
    pks = zeros(1,spwno);
    figure
    wn = round(1*sr);    % 2 sec. window (1 sec. on both sides)
    z = zeros(1,2*wn);
    for k = 1:spwno
        lfeeg = feeg(bg(k):fin(k));
        pk = find(lfeeg==max(lfeeg));
        cnt = bg(k) + pk(1);
        pks(k) = cnt;
        if cnt + wn > length(feeg)
            continue
        end
        loceeg = feeg(cnt-wn:cnt+wn);
        S1 = subplot(3,1,1);
        hold on
        plot(loceeg)
        xlim([0 2*wn])
        subplot(3,1,2)
        hold on
        locvdisc = (vdisc(vdisc>eeg_times(cnt-wn)&vdisc<eeg_times(cnt+wn))) - eeg_times(cnt) + wn / sr;
        lwn = round(locvdisc*sr);
        lwn(lwn<1) = [];
        z(lwn) = 1;
        plot(z)
        xlim([0 2*wn])
    end
    edges = 0:0.1:2;
    cnts = (edges(1:end-1) + edges(2:end)) / 2;
    z2 = histc(find(z)/sr,edges);
    z2 = z2(1:end-1);
    subplot(3,1,3)
    bar(cnts,z2)
    ts = mtz{o,1};
    ts(ts=='_') = ' ';
    axes(S1)
    title(ts)
    fns = [ts '_RIPPLEb.fig'];  % save
    saveas(gcf,fns)
    
% Ripples - centered on ripple end
    figure
    wn = round(1*sr);    % 2 sec. window (1 sec. on both sides)
    z = zeros(1,2*wn);
    for k = 1:spwno
        cnt = fin(k);
        if cnt + wn > length(feeg)
            continue
        end
        loceeg = feeg(cnt-wn:cnt+wn);
        S1 = subplot(3,1,1);
        hold on
        plot(loceeg)
        xlim([0 2*wn])
        subplot(3,1,2)
        hold on
        locvdisc = (vdisc(vdisc>eeg_times(cnt-wn)&vdisc<eeg_times(cnt+wn))) - eeg_times(cnt) + wn / sr;
        lwn = round(locvdisc*sr);
        lwn(lwn<1) = [];
        z(lwn) = 1;
        plot(z)
        xlim([0 2*wn])
    end
    edges = 0:0.1:2;
    cnts = (edges(1:end-1) + edges(2:end)) / 2;
    z2 = histc(find(z)/sr,edges);
    z2 = z2(1:end-1);
    subplot(3,1,3)
    bar(cnts,z2)
    ts = mtz{o,1};
    ts(ts=='_') = ' ';
    axes(S1)
    title(ts)
    fns = [ts '_RIPPLEc.fig'];  % save
    saveas(gcf,fns)
    
% Ripples - centered on ripple start
    figure
    wn = round(1*sr);    % 2 sec. window (1 sec. on both sides)
    z = zeros(1,2*wn);
    for k = 1:spwno
        cnt = bg(k);
        if cnt + wn > length(feeg)
            continue
        end
        loceeg = feeg(cnt-wn:cnt+wn);
        S1 = subplot(3,1,1);
        hold on
        plot(loceeg)
        xlim([0 2*wn])
        subplot(3,1,2)
        hold on
        locvdisc = (vdisc(vdisc>eeg_times(cnt-wn)&vdisc<eeg_times(cnt+wn))) - eeg_times(cnt) + wn / sr;
        lwn = round(locvdisc*sr);
        lwn(lwn<1) = [];
        z(lwn) = 1;
        plot(z)
        xlim([0 2*wn])
    end
    edges = 0:0.1:2;
    cnts = (edges(1:end-1) + edges(2:end)) / 2;
    z2 = histc(find(z)/sr,edges);
    z2 = z2(1:end-1);
    subplot(3,1,3)
    bar(cnts,z2)
    ts = mtz{o,1};
    ts(ts=='_') = ' ';
    axes(S1)
    title(ts)
    fns = [ts '_RIPPLEa.fig'];  % save
    saveas(gcf,fns)
    
% Normalized ripple histogram
    term = spwno;   % terminate if spiking ceased
    pterm = 1;      % begin when spiking starts
    for k = 1:spwno
        be = linterp1q((1:length(eeg_times))',eeg_times',bg(k));
        vdpre = vdisc(find(vdisc<be,1,'first'));
        if isempty(vdpre)
            pterm = k + 1;
        end
        fe = linterp1q((1:length(eeg_times))',eeg_times',fin(k));
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
        nm = zeros(1,length(edges)-1);
        for t = 1:length(nm)
            x1 = edges(t);
            x2 = edges(t+1);
            yy = linterp1q((1:length(eeg_times))',eeg_times',[x1; x2]);
            nm(t) = length(find(vdisc>=yy(1)&vdisc<yy(2)));
        end
        intedge = linterp1q((1:length(eeg_times))',eeg_times',edges');
        [nm2,xout] = histc(vdisc,intedge);
        isequal(nm,nm2(1:end-1)')
        nmn(k,:) = nm ./ [ones(1,12)*ivm ones(1,4)*iv1 ones(1,4)*iv2 ones(1,12)*ivm] * sr; % in spikes/sec
        
        vdpre = vdisc(find(vdisc<intedge(1),1,'last'));    % permutation test
        vdpost = vdisc(find(vdisc>intedge(end),1,'first'));
        lvd = vdisc(vdisc>=vdpre&vdisc<=vdpost);
        dlvd = diff(lvd);
        lisi = length(dlvd);
        for kk = 1:nrp
            rp = randperm(lisi);
            vdt = vdpre + cumsum(dlvd(rp));
            rnm = histc(vdt,intedge)';
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
    title(ts)
    fns = [ts '_RIPPLEk.fig'];  % save (k after Klausberger)
    saveas(gcf,fns)
    
    close all
end
cd(mm)

% -------------------------------------------------------------------------
function yi = linterp1q(x,y,xi)

yi = zeros(1,length(xi));
for k = 1:length(xi)    % linear interpolation
    inx1 = find(x<xi(k),1,'last');
    inx2 = inx1 + 1;
    yi(k) = y(inx1) + (y(inx2) - y(inx1)) * (xi(k)-x(inx1)) / (x(inx2) - x(inx1));
end
    