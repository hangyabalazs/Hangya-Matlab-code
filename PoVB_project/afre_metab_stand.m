function afre_metab_stand(inpdir1)
%AFRE_METAB_STAND   Analysis of slow and fast cycles in restricted frequency bands.
%   AFRE_METAB_STAND(DR) calculates and saves phase histograms for slow 
%   (1.1-1.6 Hz) and fast (2-2.5 Hz) cycles in two cases: (i) analysis is
%   restricted to those segments where EEG frequency is between 1.1 and 1.6
%   Hz on average, (ii) analysis is restricted to those segments where EEG
%   frequency is between 2 and 2.5 Hz on average. Percentage of missed
%   cycles and cycle length distributions are also saved. EEG is
%   standardized before phase calculations. Input directory should be given
%   as an argument (DR).
%
%   See also AFRE_RESTRICT_STAND.

% Input argument check
error(nargchk(1,1,nargin))

% Directories
global DATAPATH
inpdir_bas = [inpdir1 'bas\'];
inpdir_bic = [inpdir1 'bic\'];
inpdir2 = [DATAPATH 'Andi\Ketxyl\Cluster\mat2\'];   % burst analysis data
resdir1 = [DATAPATH 'Andi\Ketxyl\FreBandRestrict_metab\'];
mm = pwd;
dbstop if error

% Filelist
[files1_bas files_short1_bas] = filelist(inpdir_bas);
[files1_bic files_short1_bic] = filelist(inpdir_bic);
[files2 files_short2] = filelist2(inpdir2);
files_short_bas = intersect(files_short1_bas,files_short2);
files_short_bic = intersect(files_short1_bic,files_short2);
sf_bas = length(files_short_bas);
sf_bic = length(files_short_bic);

% Progress indicator
[wb,awb1,awb2] = waitbar2([0 0],'Running AFRE METAB...');
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
aang_fs1 = [];
aang_fs_slow1 = [];
aang_fs_fast1 = [];
aang_sp1 = [];
aang_sp_slow1 = [];
aang_sp_fast1 = [];
aang_as1 = [];
aang_as_slow1 = [];
aang_as_fast1 = [];
aang_afsp1 = [];
aang_afsp_slow1 = [];
aang_afsp_fast1 = [];
aang_fs2 = [];
aang_fs_slow2 = [];
aang_fs_fast2 = [];
aang_sp2 = [];
aang_sp_slow2 = [];
aang_sp_fast2 = [];
aang_as2 = [];
aang_as_slow2 = [];
aang_as_fast2 = [];
aang_afsp2 = [];
aang_afsp_slow2 = [];
aang_afsp_fast2 = [];
cycinfo_fs1 = struct('all',[],'all_slow',[],'all_fast',[],'missed',[],'missed_slow',[],'missed_fast',[]);
cycinfo_as1 = struct('all',[],'all_slow',[],'all_fast',[],'missed',[],'missed_slow',[],'missed_fast',[]);
cycinfo_sp1 = struct('all',[],'all_slow',[],'all_fast',[],'missed',[],'missed_slow',[],'missed_fast',[]);
cycinfo_afsp1 = struct('all',[],'all_slow',[],'all_fast',[],'missed',[],'missed_slow',[],'missed_fast',[]);
cycinfo_fs2 = struct('all',[],'all_slow',[],'all_fast',[],'missed',[],'missed_slow',[],'missed_fast',[]);
cycinfo_as2 = struct('all',[],'all_slow',[],'all_fast',[],'missed',[],'missed_slow',[],'missed_fast',[]);
cycinfo_sp2 = struct('all',[],'all_slow',[],'all_fast',[],'missed',[],'missed_slow',[],'missed_fast',[]);
cycinfo_afsp2 = struct('all',[],'all_slow',[],'all_fast',[],'missed',[],'missed_slow',[],'missed_fast',[]);
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
        freq = 1 / cyclen * 1000;
        if 1.1 < freq & freq < 1.6        % restrict frequency band
            lvb1 = vb1(vb1>ind1(k)&vb1<ind2(k)) - ind1(k);
            lvb1 = round(lvb1/const);    % burst first spikes, downsample unit on 1000 Hz
            [paang_fs paang_fs_slow paang_fs_fast cycinfo_fs] = laphase_stand(lfeeg,lahee,lvb1,dsr);    % PHASE - burst first spikes
            aang_fs1 = [aang_fs1 paang_fs];
            aang_fs_slow1 = [aang_fs_slow1 paang_fs_slow];
            aang_fs_fast1 = [aang_fs_fast1 paang_fs_fast];
            cycinfo_fs1.all = [cycinfo_fs1.all cycinfo_fs.all];
            cycinfo_fs1.all_slow = [cycinfo_fs1.all_slow cycinfo_fs.all_slow];
            cycinfo_fs1.all_fast = [cycinfo_fs1.all_fast cycinfo_fs.all_fast];
            cycinfo_fs1.missed = [cycinfo_fs1.missed cycinfo_fs.missed];
            cycinfo_fs1.missed_slow = [cycinfo_fs1.missed_slow cycinfo_fs.missed_slow];
            cycinfo_fs1.missed_fast = [cycinfo_fs1.missed_fast cycinfo_fs.missed_fast];
            
            lsspo = sspo(sspo>ind1(k)&sspo<ind2(k)) - ind1(k);
            lsspo = round(lsspo/const);    % single spikes, downsample unit on 1000 Hz
            [paang_sp paang_sp_slow paang_sp_fast cycinfo_sp] = laphase_stand(lfeeg,lahee,lsspo,dsr);    % PHASE - single spikes
            aang_sp1 = [aang_sp1 paang_sp];
            aang_sp_slow1 = [aang_sp_slow1 paang_sp_slow];
            aang_sp_fast1 = [aang_sp_fast1 paang_sp_fast];
            cycinfo_sp1.all = [cycinfo_sp1.all cycinfo_sp.all];
            cycinfo_sp1.all_slow = [cycinfo_sp1.all_slow cycinfo_sp.all_slow];
            cycinfo_sp1.all_fast = [cycinfo_sp1.all_fast cycinfo_sp.all_fast];
            cycinfo_sp1.missed = [cycinfo_sp1.missed cycinfo_sp.missed];
            cycinfo_sp1.missed_slow = [cycinfo_sp1.missed_slow cycinfo_sp.missed_slow];
            cycinfo_sp1.missed_fast = [cycinfo_sp1.missed_fast cycinfo_sp.missed_fast];
                        
            [paang_as paang_as_slow paang_as_fast cycinfo_as] = laphase_stand(lfeeg,lahee,vdisc2,dsr);    % PHASE - all spikes
            aang_as1 = [aang_as1 paang_as];
            aang_as_slow1 = [aang_as_slow1 paang_as_slow];
            aang_as_fast1 = [aang_as_fast1 paang_as_fast];
            cycinfo_as1.all = [cycinfo_as1.all cycinfo_as.all];
            cycinfo_as1.all_slow = [cycinfo_as1.all_slow cycinfo_as.all_slow];
            cycinfo_as1.all_fast = [cycinfo_as1.all_fast cycinfo_as.all_fast];
            cycinfo_as1.missed = [cycinfo_as1.missed cycinfo_as.missed];
            cycinfo_as1.missed_slow = [cycinfo_as1.missed_slow cycinfo_as.missed_slow];
            cycinfo_as1.missed_fast = [cycinfo_as1.missed_fast cycinfo_as.missed_fast];
            
            lafsp = afsp(afsp>ind1(k)&afsp<ind2(k)) - ind1(k);
            lafsp = round(lafsp/const);    % all first spikes, downsample unit on 1000 Hz
            [paang_afsp paang_afsp_slow paang_afsp_fast cycinfo_afsp] = laphase_stand(lfeeg,lahee,lafsp,dsr);    % PHASE - all first spikes
            aang_afsp1 = [aang_afsp1 paang_afsp];
            aang_afsp_slow1 = [aang_afsp_slow1 paang_afsp_slow];
            aang_afsp_fast1 = [aang_afsp_fast1 paang_afsp_fast];
            cycinfo_afsp1.all = [cycinfo_afsp1.all cycinfo_afsp.all];
            cycinfo_afsp1.all_slow = [cycinfo_afsp1.all_slow cycinfo_afsp.all_slow];
            cycinfo_afsp1.all_fast = [cycinfo_afsp1.all_fast cycinfo_afsp.all_fast];
            cycinfo_afsp1.missed = [cycinfo_afsp1.missed cycinfo_afsp.missed];
            cycinfo_afsp1.missed_slow = [cycinfo_afsp1.missed_slow cycinfo_afsp.missed_slow];
            cycinfo_afsp1.missed_fast = [cycinfo_afsp1.missed_fast cycinfo_afsp.missed_fast];
        end
        if 2 < freq & freq < 2.5        % restrict frequency band
            lvb1 = vb1(vb1>ind1(k)&vb1<ind2(k)) - ind1(k);
            lvb1 = round(lvb1/const);    % burst first spikes, downsample unit on 1000 Hz
            [paang_fs paang_fs_slow paang_fs_fast cycinfo_fs] = laphase_stand(lfeeg,lahee,lvb1,dsr);    % PHASE - burst first spikes
            aang_fs2 = [aang_fs2 paang_fs];
            aang_fs_slow2 = [aang_fs_slow2 paang_fs_slow];
            aang_fs_fast2 = [aang_fs_fast2 paang_fs_fast];
            cycinfo_fs2.all = [cycinfo_fs2.all cycinfo_fs.all];
            cycinfo_fs2.all_slow = [cycinfo_fs2.all_slow cycinfo_fs.all_slow];
            cycinfo_fs2.all_fast = [cycinfo_fs2.all_fast cycinfo_fs.all_fast];
            cycinfo_fs2.missed = [cycinfo_fs2.missed cycinfo_fs.missed];
            cycinfo_fs2.missed_slow = [cycinfo_fs2.missed_slow cycinfo_fs.missed_slow];
            cycinfo_fs2.missed_fast = [cycinfo_fs2.missed_fast cycinfo_fs.missed_fast];
            
            lsspo = sspo(sspo>ind1(k)&sspo<ind2(k)) - ind1(k);
            lsspo = round(lsspo/const);    % single spikes, downsample unit on 1000 Hz
            [paang_sp paang_sp_slow paang_sp_fast cycinfo_sp] = laphase_stand(lfeeg,lahee,lsspo,dsr);    % PHASE - single spikes
            aang_sp2 = [aang_sp2 paang_sp];
            aang_sp_slow2 = [aang_sp_slow2 paang_sp_slow];
            aang_sp_fast2 = [aang_sp_fast2 paang_sp_fast];
            cycinfo_sp2.all = [cycinfo_sp2.all cycinfo_sp.all];
            cycinfo_sp2.all_slow = [cycinfo_sp2.all_slow cycinfo_sp.all_slow];
            cycinfo_sp2.all_fast = [cycinfo_sp2.all_fast cycinfo_sp.all_fast];
            cycinfo_sp2.missed = [cycinfo_sp2.missed cycinfo_sp.missed];
            cycinfo_sp2.missed_slow = [cycinfo_sp2.missed_slow cycinfo_sp.missed_slow];
            cycinfo_sp2.missed_fast = [cycinfo_sp2.missed_fast cycinfo_sp.missed_fast];
                        
            [paang_as paang_as_slow paang_as_fast cycinfo_as] = laphase_stand(lfeeg,lahee,vdisc2,dsr);    % PHASE - all spikes
            aang_as2 = [aang_as2 paang_as];
            aang_as_slow2 = [aang_as_slow2 paang_as_slow];
            aang_as_fast2 = [aang_as_fast2 paang_as_fast];
            cycinfo_as2.all = [cycinfo_as2.all cycinfo_as.all];
            cycinfo_as2.all_slow = [cycinfo_as2.all_slow cycinfo_as.all_slow];
            cycinfo_as2.all_fast = [cycinfo_as2.all_fast cycinfo_as.all_fast];
            cycinfo_as2.missed = [cycinfo_as2.missed cycinfo_as.missed];
            cycinfo_as2.missed_slow = [cycinfo_as2.missed_slow cycinfo_as.missed_slow];
            cycinfo_as2.missed_fast = [cycinfo_as2.missed_fast cycinfo_as.missed_fast];
            
            lafsp = afsp(afsp>ind1(k)&afsp<ind2(k)) - ind1(k);
            lafsp = round(lafsp/const);    % all first spikes, downsample unit on 1000 Hz
            [paang_afsp paang_afsp_slow paang_afsp_fast cycinfo_afsp] = laphase_stand(lfeeg,lahee,lafsp,dsr);    % PHASE - all first spikes
            aang_afsp2 = [aang_afsp2 paang_afsp];
            aang_afsp_slow2 = [aang_afsp_slow2 paang_afsp_slow];
            aang_afsp_fast2 = [aang_afsp_fast2 paang_afsp_fast];
            cycinfo_afsp2.all = [cycinfo_afsp2.all cycinfo_afsp.all];
            cycinfo_afsp2.all_slow = [cycinfo_afsp2.all_slow cycinfo_afsp.all_slow];
            cycinfo_afsp2.all_fast = [cycinfo_afsp2.all_fast cycinfo_afsp.all_fast];
            cycinfo_afsp2.missed = [cycinfo_afsp2.missed cycinfo_afsp.missed];
            cycinfo_afsp2.missed_slow = [cycinfo_afsp2.missed_slow cycinfo_afsp.missed_slow];
            cycinfo_afsp2.missed_fast = [cycinfo_afsp2.missed_fast cycinfo_afsp.missed_fast];
        end
        waitbar2([(o-1)/sf k/lenr],wb,strw);
    end
end
if ~isempty(aang_afsp1)     % lower frequency band
    [aang_fs1 ang_fs1 mvl_fs1 ftm_fs1 n_fs1 nm_fs1] = phist(aang_fs1);
    [aang_fs_slow1 ang_fs_slow1 mvl_fs_slow1 ftm_fs_slow1 n_fs_slow1 nm_fs_slow1] = phist(aang_fs_slow1);
    [aang_fs_fast1 ang_fs_fast1 mvl_fs_fast1 ftm_fs_fast1 n_fs_fast1 nm_fs_fast1] = phist(aang_fs_fast1);
    
    [aang_sp1 ang_sp1 mvl_sp1 ftm_sp1 n_sp1 nm_sp1] = phist(aang_sp1);
    [aang_sp_slow1 ang_sp_slow1 mvl_sp_slow1 ftm_sp_slow1 n_sp_slow1 nm_sp_slow1] = phist(aang_sp_slow1);
    [aang_sp_fast1 ang_sp_fast1 mvl_sp_fast1 ftm_sp_fast1 n_sp_fast1 nm_sp_fast1] = phist(aang_sp_fast1);
    
    [aang_as1 ang_as1 mvl_as1 ftm_as1 n_as1 nm_as1] = phist(aang_as1);
    [aang_as_slow1 ang_as_slow1 mvl_as_slow1 ftm_as_slow1 n_as_slow1 nm_as_slow1] = phist(aang_as_slow1);
    [aang_as_fast1 ang_as_fast1 mvl_as_fast1 ftm_as_fast1 n_as_fast1 nm_as_fast1] = phist(aang_as_fast1);
    
    [aang_afsp1 ang_afsp1 mvl_afsp1 ftm_afsp1 n_afsp1 nm_afsp1] = phist(aang_afsp1);
    [aang_afsp_slow1 ang_afsp_slow1 mvl_afsp_slow1 ftm_afsp_slow1 n_afsp_slow1 nm_afsp_slow1] = phist(aang_afsp_slow1);
    [aang_afsp_fast1 ang_afsp_fast1 mvl_afsp_fast1 ftm_afsp_fast1 n_afsp_fast1 nm_afsp_fast1] = phist(aang_afsp_fast1);
    
    cd(resdir1)
    H = figure;      % phase histograms
    set(gcf,'Position',[1 31 1278 920])
    subplot(2,2,1)      % all spikes
    bar(cnts,nm_as1'/length(aang_as1))
    cmps = strread(fname,'%s','delimiter','_');
    titlestr = [];
    for tt = 1:length(cmps)
        titlestr = [titlestr ' ' cmps{tt}];
    end
    title(gca,[titlestr ' slow segments, all spikes'])
    y_lim = ylim;
    axis([-200 200 y_lim(1) y_lim(2)])
    str = ['\it{Mean angle: }' '\bf ' num2str(ang_as1)];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','red')
    str = ['\it{Mean resultant length: }' '\bf ' num2str(mvl_as1)];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/12,str,'Color','red')
    str = ['\it{n: }' '\bf ' num2str(n_as1)];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/6,str,'Color','red')
    subplot(2,2,2)      % all spikes - slow cycles
    bar(cnts,nm_as_slow1'/length(aang_as_slow1))
    title(gca,[titlestr ' slow segments, all spikes of slow cycles'])
    y_lim = ylim;
    axis([-200 200 y_lim(1) y_lim(2)])
    str = ['\it{Mean angle: }' '\bf ' num2str(ang_as_slow1)];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','red')
    str = ['\it{Mean resultant length: }' '\bf ' num2str(mvl_as_slow1)];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/12,str,'Color','red')
    str = ['\it{n: }' '\bf ' num2str(n_as_slow1)];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/6,str,'Color','red')
    subplot(2,2,3)      % all spikes - fast cycles
    bar(cnts,nm_as_fast1'/length(aang_as_fast1))
    title(gca,[titlestr ' slow segments, all spikes of fast cycles'])
    y_lim = ylim;
    axis([-200 200 y_lim(1) y_lim(2)])
    str = ['\it{Mean angle: }' '\bf ' num2str(ang_as_fast1)];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','red')
    str = ['\it{Mean resultant length: }' '\bf ' num2str(mvl_as_fast1)];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/12,str,'Color','red')
    str = ['\it{n: }' '\bf ' num2str(n_as_fast1)];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/6,str,'Color','red')
    str = {['% of missed cycles: ' num2str(length(cycinfo_as1.missed)/length(cycinfo_as1.all))]...
        ['% of missed slow cycles: ' num2str(length(cycinfo_as1.missed_slow)/length(cycinfo_as1.all_slow))]...
        ['% of missed fast cycles: ' num2str(length(cycinfo_as1.missed_fast)/length(cycinfo_as1.all_fast))]};
    uicontrol('Style','text','Unit','normalized','Position',...
        [0.56 0.10 0.35 0.35],'FontSize',12,'HorizontalAlignment',...
        'left','String',str,'BackgroundColor',[1 1 1]);
    dbclear if error
    fns = [fname(1:end-4) '_AS1.fig'];
    saveas(H,fns)
    fns = [fname(1:end-4) '_AS1.jpg'];
    saveas(H,fns)
    
    H = figure;      % phase histograms
    set(gcf,'Position',[1 31 1278 920])
    subplot(2,2,1)      % all spikes
    bar(cnts,nm_sp1'/length(aang_sp1))
    title(gca,[titlestr ' slow segments, single spikes'])
    y_lim = ylim;
    axis([-200 200 y_lim(1) y_lim(2)])
    str = ['\it{Mean angle: }' '\bf ' num2str(ang_sp1)];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','red')
    str = ['\it{Mean resultant length: }' '\bf ' num2str(mvl_sp1)];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/12,str,'Color','red')
    str = ['\it{n: }' '\bf ' num2str(n_sp1)];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/6,str,'Color','red')
    subplot(2,2,2)      % all spikes - slow cycles
    bar(cnts,nm_sp_slow1'/length(aang_sp_slow1))
    title(gca,[titlestr ' slow segments, single spikes of slow cycles'])
    y_lim = ylim;
    axis([-200 200 y_lim(1) y_lim(2)])
    str = ['\it{Mean angle: }' '\bf ' num2str(ang_sp_slow1)];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','red')
    str = ['\it{Mean resultant length: }' '\bf ' num2str(mvl_sp_slow1)];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/12,str,'Color','red')
    str = ['\it{n: }' '\bf ' num2str(n_sp_slow1)];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/6,str,'Color','red')
    subplot(2,2,3)      % all spikes - fast cycles
    bar(cnts,nm_sp_fast1'/length(aang_sp_fast1))
    title(gca,[titlestr ' slow segments, single spikes of fast cycles'])
    y_lim = ylim;
    axis([-200 200 y_lim(1) y_lim(2)])
    str = ['\it{Mean angle: }' '\bf ' num2str(ang_sp_fast1)];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','red')
    str = ['\it{Mean resultant length: }' '\bf ' num2str(mvl_sp_fast1)];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/12,str,'Color','red')
    str = ['\it{n: }' '\bf ' num2str(n_sp_fast1)];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/6,str,'Color','red')
    str = {['% of missed cycles: ' num2str(length(cycinfo_sp1.missed)/length(cycinfo_sp1.all))]...
        ['% of missed slow cycles: ' num2str(length(cycinfo_sp1.missed_slow)/length(cycinfo_sp1.all_slow))]...
        ['% of missed fast cycles: ' num2str(length(cycinfo_sp1.missed_fast)/length(cycinfo_sp1.all_fast))]};
    uicontrol('Style','text','Unit','normalized','Position',...
        [0.56 0.10 0.35 0.35],'FontSize',12,'HorizontalAlignment',...
        'left','String',str,'BackgroundColor',[1 1 1]);
    fns = [fname(1:end-4) '_SP1.fig'];
    saveas(H,fns)
    fns = [fname(1:end-4) '_SP1.jpg'];
    saveas(H,fns)
    
    H = figure;      % phase histograms
    set(gcf,'Position',[1 31 1278 920])
    subplot(2,2,1)      % all spikes
    bar(cnts,nm_fs1'/length(aang_fs1))
    title(gca,[titlestr ' slow segments, burst first spikes'])
    y_lim = ylim;
    axis([-200 200 y_lim(1) y_lim(2)])
    str = ['\it{Mean angle: }' '\bf ' num2str(ang_fs1)];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','red')
    str = ['\it{Mean resultant length: }' '\bf ' num2str(mvl_fs1)];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/12,str,'Color','red')
    str = ['\it{n: }' '\bf ' num2str(n_fs1)];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/6,str,'Color','red')
    subplot(2,2,2)      % all spikes - slow cycles
    bar(cnts,nm_fs_slow1'/length(aang_fs_slow1))
    title(gca,[titlestr ' slow segments, burst first spikes of slow cycles'])
    y_lim = ylim;
    axis([-200 200 y_lim(1) y_lim(2)])
    str = ['\it{Mean angle: }' '\bf ' num2str(ang_fs_slow1)];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','red')
    str = ['\it{Mean resultant length: }' '\bf ' num2str(mvl_fs_slow1)];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/12,str,'Color','red')
    str = ['\it{n: }' '\bf ' num2str(n_fs_slow1)];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/6,str,'Color','red')
    subplot(2,2,3)      % all spikes - fast cycles
    bar(cnts,nm_fs_fast1'/length(aang_fs_fast1))
    title(gca,[titlestr ' slow segments, burst first spikes of fast cycles'])
    y_lim = ylim;
    axis([-200 200 y_lim(1) y_lim(2)])
    str = ['\it{Mean angle: }' '\bf ' num2str(ang_fs_fast1)];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','red')
    str = ['\it{Mean resultant length: }' '\bf ' num2str(mvl_fs_fast1)];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/12,str,'Color','red')
    str = ['\it{n: }' '\bf ' num2str(n_fs_fast1)];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/6,str,'Color','red')
    str = {['% of missed cycles: ' num2str(length(cycinfo_fs1.missed)/length(cycinfo_fs1.all))]...
        ['% of missed slow cycles: ' num2str(length(cycinfo_fs1.missed_slow)/length(cycinfo_fs1.all_slow))]...
        ['% of missed fast cycles: ' num2str(length(cycinfo_fs1.missed_fast)/length(cycinfo_fs1.all_fast))]};
    uicontrol('Style','text','Unit','normalized','Position',...
        [0.56 0.10 0.35 0.35],'FontSize',12,'HorizontalAlignment',...
        'left','String',str,'BackgroundColor',[1 1 1]);
    fns = [fname(1:end-4) '_FS1.fig'];
    saveas(H,fns)
    fns = [fname(1:end-4) '_FS1.jpg'];
    saveas(H,fns)
    
    H = figure;      % phase histograms
    set(gcf,'Position',[1 31 1278 920])
    subplot(2,2,1)      % all spikes
    bar(cnts,nm_afsp1'/length(aang_afsp1))
    title(gca,[titlestr ' slow segments, all first spikes'])
    y_lim = ylim;
    axis([-200 200 y_lim(1) y_lim(2)])
    str = ['\it{Mean angle: }' '\bf ' num2str(ang_afsp1)];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','red')
    str = ['\it{Mean resultant length: }' '\bf ' num2str(mvl_afsp1)];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/12,str,'Color','red')
    str = ['\it{n: }' '\bf ' num2str(n_afsp1)];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/6,str,'Color','red')
    subplot(2,2,2)      % all spikes - slow cycles
    bar(cnts,nm_afsp_slow1'/length(aang_afsp_slow1))
    title(gca,[titlestr ' slow segments, all first spikes of slow cycles'])
    y_lim = ylim;
    axis([-200 200 y_lim(1) y_lim(2)])
    str = ['\it{Mean angle: }' '\bf ' num2str(ang_afsp_slow1)];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','red')
    str = ['\it{Mean resultant length: }' '\bf ' num2str(mvl_afsp_slow1)];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/12,str,'Color','red')
    str = ['\it{n: }' '\bf ' num2str(n_afsp_slow1)];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/6,str,'Color','red')
    subplot(2,2,3)      % all spikes - fast cycles
    bar(cnts,nm_afsp_fast1'/length(aang_afsp_fast1))
    title(gca,[titlestr ' slow segments, all first spikes of fast cycles'])
    y_lim = ylim;
    axis([-200 200 y_lim(1) y_lim(2)])
    str = ['\it{Mean angle: }' '\bf ' num2str(ang_afsp_fast1)];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','red')
    str = ['\it{Mean resultant length: }' '\bf ' num2str(mvl_afsp_fast1)];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/12,str,'Color','red')
    str = ['\it{n: }' '\bf ' num2str(n_afsp_fast1)];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/6,str,'Color','red')
    str = {['% of missed cycles: ' num2str(length(cycinfo_afsp1.missed)/length(cycinfo_afsp1.all))]...
        ['% of missed slow cycles: ' num2str(length(cycinfo_afsp1.missed_slow)/length(cycinfo_afsp1.all_slow))]...
        ['% of missed fast cycles: ' num2str(length(cycinfo_afsp1.missed_fast)/length(cycinfo_afsp1.all_fast))]};
    uicontrol('Style','text','Unit','normalized','Position',...
        [0.56 0.10 0.35 0.35],'FontSize',12,'HorizontalAlignment',...
        'left','String',str,'BackgroundColor',[1 1 1]);
    fns = [fname(1:end-4) '_AFSP1.fig'];
    saveas(H,fns)
    fns = [fname(1:end-4) '_AFSP1.jpg'];
    saveas(H,fns)
    
    cbox(cycinfo_afsp1.all,cycinfo_afsp1.missed)
    fns = [fname(1:end-4) '_CYCLEN1.fig'];
    saveas(gcf,fns)
    fns = [fname(1:end-4) '_CYCLEN1.jpg'];
    saveas(gcf,fns)
end
if ~isempty(aang_afsp2)     % higher frequency band
    [aang_fs2 ang_fs2 mvl_fs2 ftm_fs2 n_fs2 nm_fs2] = phist(aang_fs2);
    [aang_fs_slow2 ang_fs_slow2 mvl_fs_slow2 ftm_fs_slow2 n_fs_slow2 nm_fs_slow2] = phist(aang_fs_slow2);
    [aang_fs_fast2 ang_fs_fast2 mvl_fs_fast2 ftm_fs_fast2 n_fs_fast2 nm_fs_fast2] = phist(aang_fs_fast2);
    
    [aang_sp2 ang_sp2 mvl_sp2 ftm_sp2 n_sp2 nm_sp2] = phist(aang_sp2);
    [aang_sp_slow2 ang_sp_slow2 mvl_sp_slow2 ftm_sp_slow2 n_sp_slow2 nm_sp_slow2] = phist(aang_sp_slow2);
    [aang_sp_fast2 ang_sp_fast2 mvl_sp_fast2 ftm_sp_fast2 n_sp_fast2 nm_sp_fast2] = phist(aang_sp_fast2);
    
    [aang_as2 ang_as2 mvl_as2 ftm_as2 n_as2 nm_as2] = phist(aang_as2);
    [aang_as_slow2 ang_as_slow2 mvl_as_slow2 ftm_as_slow2 n_as_slow2 nm_as_slow2] = phist(aang_as_slow2);
    [aang_as_fast2 ang_as_fast2 mvl_as_fast2 ftm_as_fast2 n_as_fast2 nm_as_fast2] = phist(aang_as_fast2);
    
    [aang_afsp2 ang_afsp2 mvl_afsp2 ftm_afsp2 n_afsp2 nm_afsp2] = phist(aang_afsp2);
    [aang_afsp_slow2 ang_afsp_slow2 mvl_afsp_slow2 ftm_afsp_slow2 n_afsp_slow2 nm_afsp_slow2] = phist(aang_afsp_slow2);
    [aang_afsp_fast2 ang_afsp_fast2 mvl_afsp_fast2 ftm_afsp_fast2 n_afsp_fast2 nm_afsp_fast2] = phist(aang_afsp_fast2);
    
    cd(resdir1)
    H = figure;      % phase histograms
    set(gcf,'Position',[1 31 1278 920])
    subplot(2,2,1)      % all spikes
    bar(cnts,nm_as2'/length(aang_as2))
    cmps = strread(fname,'%s','delimiter','_');
    titlestr = [];
    for tt = 1:length(cmps)
        titlestr = [titlestr ' ' cmps{tt}];
    end
    title(gca,[titlestr ' fast segments, all spikes'])
    y_lim = ylim;
    axis([-200 200 y_lim(1) y_lim(2)])
    str = ['\it{Mean angle: }' '\bf ' num2str(ang_as2)];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','red')
    str = ['\it{Mean resultant length: }' '\bf ' num2str(mvl_as2)];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/12,str,'Color','red')
    str = ['\it{n: }' '\bf ' num2str(n_as2)];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/6,str,'Color','red')
    subplot(2,2,2)      % all spikes - slow cycles
    bar(cnts,nm_as_slow2'/length(aang_as_slow2))
    title(gca,[titlestr ' fast segments, all spikes of slow cycles'])
    y_lim = ylim;
    axis([-200 200 y_lim(1) y_lim(2)])
    str = ['\it{Mean angle: }' '\bf ' num2str(ang_as_slow2)];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','red')
    str = ['\it{Mean resultant length: }' '\bf ' num2str(mvl_as_slow2)];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/12,str,'Color','red')
    str = ['\it{n: }' '\bf ' num2str(n_as_slow2)];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/6,str,'Color','red')
    subplot(2,2,3)      % all spikes - fast cycles
    bar(cnts,nm_as_fast2'/length(aang_as_fast2))
    title(gca,[titlestr ' fast segments, all spikes of fast cycles'])
    y_lim = ylim;
    axis([-200 200 y_lim(1) y_lim(2)])
    str = ['\it{Mean angle: }' '\bf ' num2str(ang_as_fast2)];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','red')
    str = ['\it{Mean resultant length: }' '\bf ' num2str(mvl_as_fast2)];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/12,str,'Color','red')
    str = ['\it{n: }' '\bf ' num2str(n_as_fast2)];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/6,str,'Color','red')
    str = {['% of missed cycles: ' num2str(length(cycinfo_as2.missed)/length(cycinfo_as2.all))]...
        ['% of missed slow cycles: ' num2str(length(cycinfo_as2.missed_slow)/length(cycinfo_as2.all_slow))]...
        ['% of missed fast cycles: ' num2str(length(cycinfo_as2.missed_fast)/length(cycinfo_as2.all_fast))]};
    uicontrol('Style','text','Unit','normalized','Position',...
        [0.56 0.10 0.35 0.35],'FontSize',12,'HorizontalAlignment',...
        'left','String',str,'BackgroundColor',[1 1 1]);
    dbclear if error
    fns = [fname(1:end-4) '_AS2.fig'];
    saveas(H,fns)
    fns = [fname(1:end-4) '_AS2.jpg'];
    saveas(H,fns)
    
    H = figure;      % phase histograms
    set(gcf,'Position',[1 31 1278 920])
    subplot(2,2,1)      % all spikes
    bar(cnts,nm_sp2'/length(aang_sp2))
    title(gca,[titlestr ' fast segments, single spikes'])
    y_lim = ylim;
    axis([-200 200 y_lim(1) y_lim(2)])
    str = ['\it{Mean angle: }' '\bf ' num2str(ang_sp2)];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','red')
    str = ['\it{Mean resultant length: }' '\bf ' num2str(mvl_sp2)];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/12,str,'Color','red')
    str = ['\it{n: }' '\bf ' num2str(n_sp2)];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/6,str,'Color','red')
    subplot(2,2,2)      % all spikes - slow cycles
    bar(cnts,nm_sp_slow2'/length(aang_sp_slow2))
    title(gca,[titlestr ' fast segments, single spikes of slow cycles'])
    y_lim = ylim;
    axis([-200 200 y_lim(1) y_lim(2)])
    str = ['\it{Mean angle: }' '\bf ' num2str(ang_sp_slow2)];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','red')
    str = ['\it{Mean resultant length: }' '\bf ' num2str(mvl_sp_slow2)];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/12,str,'Color','red')
    str = ['\it{n: }' '\bf ' num2str(n_sp_slow2)];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/6,str,'Color','red')
    subplot(2,2,3)      % all spikes - fast cycles
    bar(cnts,nm_sp_fast2'/length(aang_sp_fast2))
    title(gca,[titlestr ' fast segments, single spikes of fast cycles'])
    y_lim = ylim;
    axis([-200 200 y_lim(1) y_lim(2)])
    str = ['\it{Mean angle: }' '\bf ' num2str(ang_sp_fast2)];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','red')
    str = ['\it{Mean resultant length: }' '\bf ' num2str(mvl_sp_fast2)];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/12,str,'Color','red')
    str = ['\it{n: }' '\bf ' num2str(n_sp_fast2)];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/6,str,'Color','red')
    str = {['% of missed cycles: ' num2str(length(cycinfo_sp2.missed)/length(cycinfo_sp2.all))]...
        ['% of missed slow cycles: ' num2str(length(cycinfo_sp2.missed_slow)/length(cycinfo_sp2.all_slow))]...
        ['% of missed fast cycles: ' num2str(length(cycinfo_sp2.missed_fast)/length(cycinfo_sp2.all_fast))]};
    uicontrol('Style','text','Unit','normalized','Position',...
        [0.56 0.10 0.35 0.35],'FontSize',12,'HorizontalAlignment',...
        'left','String',str,'BackgroundColor',[1 1 1]);
    fns = [fname(1:end-4) '_SP2.fig'];
    saveas(H,fns)
    fns = [fname(1:end-4) '_SP2.jpg'];
    saveas(H,fns)
    
    H = figure;      % phase histograms
    set(gcf,'Position',[1 31 1278 920])
    subplot(2,2,1)      % all spikes
    bar(cnts,nm_fs2'/length(aang_fs2))
    title(gca,[titlestr ' fast segments, burst first spikes'])
    y_lim = ylim;
    axis([-200 200 y_lim(1) y_lim(2)])
    str = ['\it{Mean angle: }' '\bf ' num2str(ang_fs2)];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','red')
    str = ['\it{Mean resultant length: }' '\bf ' num2str(mvl_fs2)];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/12,str,'Color','red')
    str = ['\it{n: }' '\bf ' num2str(n_fs2)];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/6,str,'Color','red')
    subplot(2,2,2)      % all spikes - slow cycles
    bar(cnts,nm_fs_slow2'/length(aang_fs_slow2))
    title(gca,[titlestr ' fast segments, burst first spikes of slow cycles'])
    y_lim = ylim;
    axis([-200 200 y_lim(1) y_lim(2)])
    str = ['\it{Mean angle: }' '\bf ' num2str(ang_fs_slow2)];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','red')
    str = ['\it{Mean resultant length: }' '\bf ' num2str(mvl_fs_slow2)];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/12,str,'Color','red')
    str = ['\it{n: }' '\bf ' num2str(n_fs_slow2)];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/6,str,'Color','red')
    subplot(2,2,3)      % all spikes - fast cycles
    bar(cnts,nm_fs_fast2'/length(aang_fs_fast2))
    title(gca,[titlestr ' fast segments, burst first spikes of fast cycles'])
    y_lim = ylim;
    axis([-200 200 y_lim(1) y_lim(2)])
    str = ['\it{Mean angle: }' '\bf ' num2str(ang_fs_fast2)];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','red')
    str = ['\it{Mean resultant length: }' '\bf ' num2str(mvl_fs_fast2)];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/12,str,'Color','red')
    str = ['\it{n: }' '\bf ' num2str(n_fs_fast2)];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/6,str,'Color','red')
    str = {['% of missed cycles: ' num2str(length(cycinfo_fs2.missed)/length(cycinfo_fs2.all))]...
        ['% of missed slow cycles: ' num2str(length(cycinfo_fs2.missed_slow)/length(cycinfo_fs2.all_slow))]...
        ['% of missed fast cycles: ' num2str(length(cycinfo_fs2.missed_fast)/length(cycinfo_fs2.all_fast))]};
    uicontrol('Style','text','Unit','normalized','Position',...
        [0.56 0.10 0.35 0.35],'FontSize',12,'HorizontalAlignment',...
        'left','String',str,'BackgroundColor',[1 1 1]);
    fns = [fname(1:end-4) '_FS2.fig'];
    saveas(H,fns)
    fns = [fname(1:end-4) '_FS2.jpg'];
    saveas(H,fns)
    
    H = figure;      % phase histograms
    set(gcf,'Position',[1 31 1278 920])
    subplot(2,2,1)      % all spikes
    bar(cnts,nm_afsp2'/length(aang_afsp2))
    title(gca,[titlestr ' fast segments, all first spikes'])
    y_lim = ylim;
    axis([-200 200 y_lim(1) y_lim(2)])
    str = ['\it{Mean angle: }' '\bf ' num2str(ang_afsp2)];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','red')
    str = ['\it{Mean resultant length: }' '\bf ' num2str(mvl_afsp2)];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/12,str,'Color','red')
    str = ['\it{n: }' '\bf ' num2str(n_afsp2)];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/6,str,'Color','red')
    subplot(2,2,2)      % all spikes - slow cycles
    bar(cnts,nm_afsp_slow2'/length(aang_afsp_slow2))
    title(gca,[titlestr ' fast segments, all first spikes of slow cycles'])
    y_lim = ylim;
    axis([-200 200 y_lim(1) y_lim(2)])
    str = ['\it{Mean angle: }' '\bf ' num2str(ang_afsp_slow2)];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','red')
    str = ['\it{Mean resultant length: }' '\bf ' num2str(mvl_afsp_slow2)];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/12,str,'Color','red')
    str = ['\it{n: }' '\bf ' num2str(n_afsp_slow2)];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/6,str,'Color','red')
    subplot(2,2,3)      % all spikes - fast cycles
    bar(cnts,nm_afsp_fast2'/length(aang_afsp_fast2))
    title(gca,[titlestr ' fast segments, all first spikes of fast cycles'])
    y_lim = ylim;
    axis([-200 200 y_lim(1) y_lim(2)])
    str = ['\it{Mean angle: }' '\bf ' num2str(ang_afsp_fast2)];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','red')
    str = ['\it{Mean resultant length: }' '\bf ' num2str(mvl_afsp_fast2)];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/12,str,'Color','red')
    str = ['\it{n: }' '\bf ' num2str(n_afsp_fast2)];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/6,str,'Color','red')
    str = {['% of missed cycles: ' num2str(length(cycinfo_afsp2.missed)/length(cycinfo_afsp2.all))]...
        ['% of missed slow cycles: ' num2str(length(cycinfo_afsp2.missed_slow)/length(cycinfo_afsp2.all_slow))]...
        ['% of missed fast cycles: ' num2str(length(cycinfo_afsp2.missed_fast)/length(cycinfo_afsp2.all_fast))]};
    uicontrol('Style','text','Unit','normalized','Position',...
        [0.56 0.10 0.35 0.35],'FontSize',12,'HorizontalAlignment',...
        'left','String',str,'BackgroundColor',[1 1 1]);
    fns = [fname(1:end-4) '_AFSP2.fig'];
    saveas(H,fns)
    fns = [fname(1:end-4) '_AFSP2.jpg'];
    saveas(H,fns)
    
    cbox(cycinfo_afsp2.all,cycinfo_afsp2.missed)
    fns = [fname(1:end-4) '_CYCLEN2.fig'];
    saveas(gcf,fns)
    fns = [fname(1:end-4) '_CYCLEN2.jpg'];
    saveas(gcf,fns)
end
close all
dbstop if error



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
cyclen = mean(cl6);   % cycle length in ms

% -------------------------------------------------------------------------
function [ang1 ang2 ang3 cycinfo] = laphase_stand(feeg,ahee,vdisc,sr)
%LAPHASE_STAND    Phase angles for unit relative to EEG.
%   [A1 A2 A3] = LAPHASE_STAND(FEEG,AHEE,VDISC,SR) calculates Hilbert phase 
%   angles (A1) for discriminated unit (VDISC) relative to filtered EEG 
%   (FEEG), when sampling frequency is given in SR and Hilbert-transform of
%   the EEG in AHEE. It also returns angle sets corresponding to slow (1.1-1.6)
%   and fast (2-2.5 Hz) cycles (A2). Cycles not fulfilling the following 2 criteria are
%   discarded: (i) EEG amp. higher then 2SD; (ii) min. 250 ms length.
%   Indices of discarded spikes of vdisc are returned in I.
%
%   [A1 A2 A3 CI] = LAPHASE_STAND(FEEG,AHEE,VDISC,SR) returns cycle
%   information in CI. The fields of the CI structure contains the length
%   of all cycles, slow/fast cycles and missed cycles.
%
%   See also HILBERT.

% Check SWS criteria:
% 1. discard cicles with EEG amp. lower then 2SD
% 2. discard cicles shorter then 250 ms
fn = find(-diff(ahee)>2*pi-0.1);
sd = std(feeg);
inx1 = find(vdisc<fn(1));
inx2 = find(vdisc<fn(1));
inx3 = find(vdisc<fn(1));
allcyc = [];
allcyc2 = [];
allcyc3 = [];
missedcyc = [];
missedcyc2 = [];
missedcyc3 = [];
for k = 1:length(fn)-1
    cl = (fn(k+1) - fn(k)) / sr;    % cycle length in sec.
    cf = 1 / cl;    % frequency of cycle
    seeg = feeg(fn(k):fn(k+1));
    axs = max(seeg) - min(seeg);
    sahee = ahee(fn(k):fn(k+1));
    cycvd = vdisc(vdisc>fn(k)&vdisc<fn(k+1));
    if (axs < 2 * sd)  | (fn(k+1) - fn(k) < 0.25 * sr)
        inx1 = [inx1 find(vdisc>fn(k)&vdisc<fn(k+1))];
    else
        allcyc = [allcyc cl];
        if isempty(cycvd)
            missedcyc = [missedcyc cl];
        end
        if cf > 1.1 & cf < 1.6
            allcyc2 = [allcyc2 cl];
            if isempty(cycvd)
                missedcyc2 = [missedcyc2 cl];
            end
        end
        if cf > 2 & cf < 2.5
            allcyc3 = [allcyc3 cl];
            if isempty(cycvd)
                missedcyc3 = [missedcyc3 cl];
            end
        end
    end
    if (axs < 2 * sd)  | (fn(k+1) - fn(k) < 0.25 * sr) | ~(cf > 1.1 & cf < 1.6)
        inx2 = [inx2 find(vdisc>fn(k)&vdisc<fn(k+1))];  % for slow cycles
    end
    if (axs < 2 * sd)  | (fn(k+1) - fn(k) < 0.25 * sr) | ~(cf > 2 & cf < 2.5)
        inx3 = [inx3 find(vdisc>fn(k)&vdisc<fn(k+1))];  % for fast cycles
    end
end
inx1 = [inx1 find(vdisc>fn(end))];
inx2 = [inx2 find(vdisc>fn(end))];
inx3 = [inx3 find(vdisc>fn(end))];
vdisc1 = vdisc;
vdisc1(inx1) = [];
ang1 = ahee(vdisc1);
vdisc2 = vdisc;
vdisc2(inx2) = [];
ang2 = ahee(vdisc2);
vdisc3 = vdisc;
vdisc3(inx3) = [];
ang3 = ahee(vdisc3);
cycinfo.all = allcyc;
cycinfo.all_slow = allcyc2;
cycinfo.all_fast = allcyc3;
cycinfo.missed = missedcyc;
cycinfo.missed_slow = missedcyc2;
cycinfo.missed_fast = missedcyc3;

% -------------------------------------------------------------------------
function [aang ang mvl ftm n nm] = phist(aang)

edges = -180:20:180;     % edges for phase histogram
if ~isempty(aang)
    n = length(aang);     % burst first spikes
    ftm = sum(exp(1).^(i*aang)) / n;    % first trigonometric moment
    ang = angle(ftm);   % mean angle
    mvl = abs(ftm);     % mean resultant length
    aang = aang * 180 / pi;
    ang = ang * 180 / pi;
    [nm,xout] = histc(aang,edges);   % phase histogram
    nm = nm(1:end-1);
else
    n = 0;
    ftm = NaN;
    mvl = NaN;
    aang = NaN;
    ang = NaN;
    nm = zeros(1,length(edges)-1);
end

% -------------------------------------------------------------------------
function H = cbox(m1,m2)

dbstop if error
if isempty(m1)
    m1 = 0;
end
if isempty(m2)
    m2 = 0;
end
H = figure;
boxplot([m1 m2],[zeros(size(m1)) ones(size(m2))],'labels',[{'all'} {'missed'}]);
[Wp_eu,Wh_eu] = b_ranksum2(m1,m2,'alpha',0.05);
if Wh_eu
    clr = 'red';
else
    clr = 'black';
end
y_lim = ylim;
x_lim = xlim;
tpos1 = x_lim(1) + (x_lim(2) - x_lim(1)) / 2;
tpos2 = y_lim(1) + (y_lim(2)-y_lim(1)) * 3 / 5;
text(tpos1,tpos2,num2str(Wp_eu),'Color',clr,'Horizontalalignment','center')
title('Cycle length distribution')
dbclear if error