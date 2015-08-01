function aphaserun_burst2_stand
%APHASERUN_BURST2_STAND    Phase analysis.
%   APHASERUN_BURST2_STAND calculates action potential phases relative to
%   slow LFP oscillation. It calculates the dependence of the phase
%   distribution on the burst parameters (intraburst spike no., intraburst
%   frequency, burst length. Phase histograms for all spikes, burst first
%   spikes, all first spikes and single spikes are plotted and saved, along
%   with the corresponding normalized Spike Triggered Averages. A single
%   von Mises probability density function, as well as the mixture of two
%   von Mises density functions are fitted on the phase distribution. EEG
%   wavelet and autosorrelation are also calculated, plotted and saved.
%   Edit program code to modify input and output directories!
%
%   APHASERUN_BURST2_STAND performs EEG standarization before calculation
%   Hilbert-transform (see APHASE_STAND).
%
%   See also APHASE_STAND, ASTANORM, WATSON2 and WATSONTWOFIT2.

% Directories
global DATAPATH
inpdir1 = 'X:\In_Vivo\_analysis\acsady_csoport\auj_istvan\mat_ket_xyl\';    % mat files
inpdir2 = [DATAPATH 'Andi\Ketxyl\Cluster\mat2\'];   % burst analysis data
resdir1 = [DATAPATH 'Andi\Ketxyl\AllPhase_stand\'];
mm = pwd;

% Filelist
[files1 files_short1] = filelist(inpdir1);
[files2 files_short2] = filelist2(inpdir2);
files_short = intersect(files_short1,files_short2);
sf = length(files_short);

% Main
for o = 1:sf
    H1 = figure;    % figures
    H2 = figure;
    H3 = figure;
    H4 = figure;
    H5 = figure;
    H6 = figure;
    H7 = figure;
    H8 = figure;
    
    fname = files_short{o};     % load
    ff = [inpdir1 fname];
    load(ff)
    eeg0 = data(:,2)';
    eeg = eeg0(1:20:end);    % downsample on 1000 Hz
    eeg2 = eeg0(1:100:end);    % downsample on 200 Hz
    clear data eeg0
    sr = 1000;
    ff2 = [inpdir2 fname(1:end-4) '_CLUST2.mat'];
    load(ff2)
    
    fst = Burst(1,:);
    vb1 = vdisc(fst);
    vb1 = round(vb1/20);    % burst first spikes
    
    [aang_fs dinx_fs] = aphase_stand(eeg,vb1,sr);    % PHASE - burst first spikes
    n_fs = length(aang_fs);
    ftm_fs = sum(exp(1).^(i*aang_fs)) / n_fs;    % first trigonometric moment
    ang_fs = angle(ftm_fs);   % mean angle
    mvl_fs = abs(ftm_fs);     % mean resultant length
    aang_fs = aang_fs * 180 / pi;
    ang_fs = ang_fs * 180 / pi;
    edges = [-180:20:180];
    cnts = (edges(1:end-1) + edges(2:end)) / 2;
    [nm_fs,xout_fs] = histc(aang_fs,edges);   % phase histogram
    nm_fs = nm_fs(1:end-1);
    
    sso = vdisc;       % single spikes
    for k = 1:size(Burst,2)
        sso(Burst(1,k):Burst(2,k)) = 0;
    end
    sspo = sso(sso>0);
    sspo = round(sspo/20);    % downsample unit on 1000 Hz
    
    [aang_sp dinx_sp] = aphase_stand(eeg,sspo,sr);    % PHASE - single spikes
    n_sp = length(aang_sp);
    ftm_sp = sum(exp(1).^(i*aang_sp)) / n_sp;    % first trigonometric moment
    ang_sp = angle(ftm_sp);   % mean angle
    mvl_sp = abs(ftm_sp);     % mean resultant length
    aang_sp = aang_sp * 180 / pi;
    ang_sp = ang_sp * 180 / pi;
    [nm_sp,xout_sp] = histc(aang_sp,edges);   % phase histogram
    nm_sp = nm_sp(1:end-1);
    
    [aang_as dinx_as] = aphase_stand(eeg,round(vdisc/20),sr);    % PHASE - all spikes
    n_as = length(aang_as);
    ftm_as = sum(exp(1).^(i*aang_as)) / n_as;    % first trigonometric moment
    ang_as = angle(ftm_as);   % mean angle
    mvl_as = abs(ftm_as);     % mean resultant length
    aang_as = aang_as * 180 / pi;
    ang_as = ang_as * 180 / pi;
    [nm_as,xout_as] = histc(aang_as,edges);   % phase histogram
    nm_as = nm_as(1:end-1);
    
    ssi = vdisc;       % allfirstspikes
    for k = 1:size(Burst,2)
        ssi(Burst(1,k)+1:Burst(2,k)) = 0;
    end
    afsp = ssi(ssi>0);
    afsp = round(afsp/20);    % downsample unit on 1000 Hz
    
    [aang_afsp dinx_afsp] = aphase_stand(eeg,afsp,sr);    % PHASE - all first spikes
    n_afsp = length(aang_afsp);
    ftm_afsp = sum(exp(1).^(i*aang_afsp)) / n_afsp;    % first trigonometric moment
    ang_afsp = angle(ftm_afsp);   % mean angle
    mvl_afsp = abs(ftm_afsp);     % mean resultant length
    aang_afsp = aang_afsp * 180 / pi;
    ang_afsp = ang_afsp * 180 / pi;
    [nm_afsp,xout_afsp] = histc(aang_afsp,edges);   % phase histogram
    nm_afsp = nm_afsp(1:end-1);
    
    ib = zeros(size(vdisc));       % intraburst spikes
    for k = 1:size(Burst,2)
        ib(Burst(1,k):Burst(2,k)) = 1;
        sso(Burst(1,k):Burst(2,k)) = 0;
    end
    ibs = vdisc(find(ib));
    ibs = round(ibs/20);    % downsample unit on 1000 Hz
    
    [aang_ibsang aang_sspoang aang_allang cycnb] = aphaseb(eeg,sspo,ibs,vb1,sr);
    n_ibsang = length(aang_ibsang);
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
    [nm_sspoang,xout_sspoang] = histc(aang_sspoang,edges);   % phase histogram
    nm_sspoang = nm_sspoang(1:end-1);
    n_allang = length(aang_allang);
    ftm_allang = sum(exp(1).^(i*aang_allang)) / n_allang;    % first trigonometric moment
    ang_allang = angle(ftm_allang);   % mean angle
    mvl_allang = abs(ftm_allang);     % mean resultant length
    aang_allang = aang_allang * 180 / pi;
    ang_allang = ang_allang * 180 / pi;
    [nm_allang,xout_allang] = histc(aang_allang,edges);   % phase histogram
    nm_allang = nm_allang(1:end-1);
    [nm_cycnb xout_cycnb] = hist(cycnb(cycnb>0),(1:10));     % distr. of burst no./cycle
    
    wn = 2 * sr;    % 2 sec. window
    [sta_fs sta_index1_fs sta_index2_fs nn_fs] = astanorm(vb1,eeg,wn);    % STA - burst first spikes
    [sta_sp sta_index1_sp sta_index2_sp nn_sp] = astanorm(sspo,eeg,wn);    % STA - single spikes
    [sta_as sta_index1_as sta_index2_as nn_as] = astanorm(round(vdisc/20),eeg,wn);    % STA - all spikes
    [sta_afsp sta_index1_afsp sta_index2_afsp nn_afsp] = astanorm(afsp,eeg,wn);    % STA - all first spikes
    time = linspace(-wn/sr/2,wn/sr/2,length(sta_fs));
    
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
    
    [pow,phase,f] = eeg_wavelet(eeg2);        % EEG WAVELET
%     disp('Truncated EEG!')
    
    wnd = 1;        % EEG autocorrelogram (1 sec. window)
    acr = xcorr(eeg,wnd*sr);
    acr(length(acr)/2+0.5) = [];
    acr = reshape(acr,length(acr)/200,200);
    sacr = sum(acr);
    
    figure(H1)      % burst parameters vs. phase
    set(gcf,'Position',[1 31 1278 920])     % maximize figure window
    subplot(2,2,1)
%     bar(xout,nm_fs/length(aang_fs))
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
    
    ibspno = IntraBurstSpikeNumber.all;
    ibfr = IntraBurstFrequency.all;
    bl = BurstLength.all;
    ibspno(dinx_fs) = [];
    ibfr(dinx_fs) = [];
    bl(dinx_fs) = [];
    subplot(2,2,2)
    plot(ibspno+rand(size(ibspno))*0.1,aang_fs,'.')     % intraburst spike number vs. phase
    hold on
    plot(ones(size(aang_sp))+rand(size(aang_sp))*1-1,aang_sp,'.')
    title(gca,'intraburst spike number')
    subplot(2,2,3)
    plot(ibfr,aang_fs,'.','Color','white')
    hold on
    title(gca,'intraburst frequency')     % intraburst freq. vs. phase
    clr = {'red' 'green' 'blue' 'magenta' 'yellow' 'cyan' ...
        [0 153 0]/256 [255 102 0]/256 [255 204 204]/256 [102 102 102]/256};
    mibs = max(ibspno);
    plot_handle = [];
    legstr = {};
    for k = 2:mibs
        str = ['let = find(ibspno==' num2str(k) ');'];
        eval(str)
        if ~isempty(let)
            plot_handle(end+1) = plot(ibfr(let),aang_fs(let),'.','Color',clr{k-1});
            legstr{end+1} = num2str(k);
        end
    end
    set(gca,'Color','black')
    [L,OBJ] = legend(plot_handle,legstr,'Location','NorthWestOutside');
    fo = findobj(OBJ,'Type','text');
    set(fo,'Color','white')
    subplot(2,2,4)
    plot(bl,aang_fs,'.','Color','white')
    hold on
    title(gca,'burst length')     % burst length vs. phase
    plot_handle = [];
    legstr = {};
    for k = 2:mibs
        str = ['let = find(ibspno==' num2str(k) ');'];
        eval(str)
        if ~isempty(let)
            plot_handle(end+1) = plot(bl(let),aang_fs(let),'.','Color',clr{k-1});
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
        str = ['let = find(ibspno==' num2str(k) ');'];
        eval(str)
        [nm,xout] = hist(aang_fs(let),edges);
        nmm = [nmm; nm];
        xoutt = [xoutt; xout(:)'];
        lele = length(let);
        ng = aang_fs(let) / 180 * pi;
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
    
    figure(H7)      % EEG wavelet
    set(gcf,'Position',[4 500 1220 420])
    subplot(1,2,1)
    imagesc(pow)
    b_rescaleaxis('Y',f)
    title(gca,[titlestr ' wavelet'])
    
    subplot(1,2,2)      % EEG autocorrelogram
    bar(linspace(-1000,1000,length(sacr)),sacr)
    title(gca,[titlestr ' autocorr.'])
    xlim([-1000 1000])
    
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
    fns = [fname(1:end-4) '_EEG.jpg'];
    saveas(H7,fns)
    fns = [fname(1:end-4) '_CYCLEFIRST.fig'];
    saveas(H8,fns)
    fns = [fname(1:end-4) '_CYCLEFIRST.jpg'];
    saveas(H8,fns)
    close all
end
cd(mm)

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
function [pow,phase,f] = eeg_wavelet(dat)

% Prepare for wavelet transformation
variance = std(dat) ^ 2;
dat = (dat - mean(dat)) / sqrt(variance) ;
n = length(dat);
dt = 1 / 200;
pad = 1;
dj = 0.08;    
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
[wave,period,scale,coi] = b_wavelet_new3(dat,dt,pad,dj,s0,j1,mother,param,mis);
pow = abs(wave).^2;
phase = angle(wave);

% -------------------------------------------------------------------------
function [ang_ibsang ang_sspoang ang_allang cycnb] = aphaseb(eeg,sspo,ibs,vb1,sr)
%APHASEB    Phase angles for unit relative to EEG.
%   [IBSANG SSPOANG] = APHASEB(EEG,SSPO,IBS,VB1,SR) calculates Hilbert 
%   phase angles for first intraburst and first single spike of each cycle
%   relative to the EEG, when sampling frequency is given in SR, single 
%   spikes are given in SSPO, intraburst spikes are given in IBS, burst 
%   first spikes are given in VB1. Cycles not fulfilling the following 2 
%   criteria are discarded: (i) EEG amp. higher then 2SD; (ii) min. 100 ms 
%   length (half wavelength of filter cutoff freq.). 
%
%   [IBSANG SSPOANG ALLANG] = APHASEB(EEG,SSPO,IBS,VB1,SR) returns the
%   phase of cycle first spikes as well.
%
%   [IBSANG SSPOANG ALLANG CYCNB] = APHASEB(EEG,SSPO,IBS,VB1,SR) returns the
%   number of bursts in each cycles.
%
%   See also HILBERT.

% Filtering EEG
nqf = sr / 2;
flt = fir1(4096,5/nqf,'low');      % lowpass filtering on 5 Hz
feeg = filtfilt(flt,1,eeg);
feeg = (feeg - mean(feeg)) / std(feeg);

% Hilbert transformation
ahee = angle(hilbert(feeg));

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