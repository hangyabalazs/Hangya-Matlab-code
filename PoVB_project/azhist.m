function azhist
%AZHIST    Z-shift histogram.
%   AZHIST calculates and plots 'ZMaxLoc' histogram. The histogram is saved
%   in the results directory.
%
%   See also AZSHIFTRUN5.

% Directories
global DATAPATH
inproot = 'X:\In_Vivo\_analysis\acsady_csoport\auj_istvan\mat_ket_xyl_control\';
resdir = [DATAPATH 'Andi\Ketxyl\Zshift\control\hist\'];
inpdir2 = [DATAPATH 'Andi\Ketxyl\Zshift\control\mat\'];
mm = pwd;
dr = dir(inproot);
inpadd = {};
for k = 3:length(dr)
    if dr(k).isdir
        inpadd{end+1} = dr(k).name;
    end
end

% ZMaxLoc histogram
sr = 20000;
for k = 1:length(inpadd)
    inpdir = [inproot inpadd{k} '\']
    files = b_filelist(inpdir);
    sf = length(files);
    zsh = [];
    for o = 1:sf
        cd(inpdir2)
        cmps = strread(files(o).name(1:end-4),'%s','delimiter','_');
        fn = [cmps{1} '_' cmps{2} '_' cmps{3}];
        ff = [fn '_ZSHIFT.mat'];
        load(ff)
        zsh = [zsh ZMaxLoc];
    end
    
    H = figure;     % plot and save
    [ZHist,ZHx] = hist(zsh/sr*1000,length(zsh)/10);
    bar(ZHx,ZHist)
    x_lim = xlim;
    y_lim = ylim;
    text(x_lim(1)+(x_lim(2)-x_lim(1))*0.3,y_lim(1)+(y_lim(2)-y_lim(1))*0.9,...
        ['mean: ' num2str(mean(zsh/sr*1000))])
    cd(resdir)
    eval(['saveas(H,''',inpadd{k},'_ZHIST.fig'')']);
    close(H)
    str = [inpadd{k} '_ZHIST'];
    save(str,'ZHist','ZHx')
end
cd(mm)