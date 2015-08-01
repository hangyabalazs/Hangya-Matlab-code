function aburstat

% Directories
global DATAPATH
inpdir1 = 'X:\In_Vivo\_analysis\acsady_csoport\auj_istvan\MAT\';
inpdir2 = [DATAPATH 'Andi\Disc\urethane\'];
inpdir3 = [DATAPATH 'Andi\Cluster\mat2\'];
resdir1 = [DATAPATH 'Andi\BurSTat\mat\'];
resdir2 = [DATAPATH 'Andi\BurSTat\fig\'];
resdir3 = [DATAPATH 'Andi\Length\'];
resdir4 = [DATAPATH 'Andi\BurSTat\xls\'];
mm = pwd;
cd(inpdir1)

% Filelist
[files1 files_short1] = filelist(inpdir1);
[files2 files_short2] = filelist2(inpdir2);
files_short = intersect(files_short1,files_short2);
sf = length(files_short);

% MAIN
sr = 20000;
H = figure;
nextcol = 1;
for o = 1:sf
    fname = files_short{o};
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
    clear data
    
%     efflen = (vdisc(end) - vdisc(1)) / 20000;
%     frate = length(vdisc) / efflen;
%     seglen = 50 / frate;        % adaptive segment length: should contain at least 50 spikes on average
    seglen = 30 * sr;        % 30 sec. long segments
    lenr = floor(len/seglen);
    ind1 = [1:seglen:len];
    ind2 = ind1 + seglen - 1;
    burstiness = zeros(1,lenr);
    vburst = vdisc(Burst);
    for k = 1:lenr
        vd = vdisc(find(vdisc>ind1(k)&vdisc<ind2(k)));
        lvb = locvburst(vburst,ind1(k),ind2(k));
        burstnum = size(lvb,2);
        intraburstiv = [];
        intraburstnum = zeros(1,burstnum);
        for j = 1:burstnum      % computing intraburstiv
            b = vdisc(find(vdisc>=lvb(1,j)&vdisc<=lvb(2,j)));
            db = diff(b);
            intraburstiv = [intraburstiv db];
            intraburstnum(j) = length(b);   %intraburst spike number
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
    end
    
% Save
    plot(burstiness);
    set(gca,'YLim',[0 1])
    fn = [resdir2 fname(1:end-4) '_BURSTINESS'];
    saveas(gcf,fn)
    
    plot(intraburstfreq,'.');
    set(gca,'YLim',[100 600])
    fn = [resdir2 fname(1:end-4) '_IBFR'];
    saveas(gcf,fn)
    
    plot(ibspno,'.');
    set(gca,'YLim',[0 10])
    fn = [resdir2 fname(1:end-4) '_IBSPNO'];
    saveas(gcf,fn)
    
    fn = [resdir1 fname(1:end-4) '_BURSTAT'];
    save(fn,'burstiness','intraburstfreq','ibspno','seglen')
    
    fn = [resdir3 fname(1:end-4) '_LEN'];
    save(fn,'len')
    
    colref1 = [colno2colref(nextcol) '1'];
    colref2 = [colno2colref(nextcol) '2'];
    xlswrite([resdir4 'burstiness.xls'],{fname},'sheet1',colref1)
    xlswrite([resdir4 'burstiness.xls'],burstiness','sheet1',colref2)
    xlswrite([resdir4 'intraburstfreq.xls'],{fname},'sheet1',colref1)
    xlswrite([resdir4 'intraburstfreq.xls'],intraburstfreq','sheet1',colref2)
    xlswrite([resdir4 'ibspno.xls'],{fname},'sheet1',colref1)
    xlswrite([resdir4 'ibspno.xls'],ibspno','sheet1',colref2)
    nextcol = nextcol + 1;
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