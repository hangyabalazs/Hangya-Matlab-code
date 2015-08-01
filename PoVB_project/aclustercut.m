function aclustercut
%ACLUSTERCUT    Save desired number of clusters and number of missed spikes.
%   ACLUSTERCUT loads output files of ACLUST, asks for desired number of
%   clusters and saves result. It plots histogram of afterburst intervals, asks
%   for number for missed spikes and saves result. Edit the program code to
%   modify the input and result directories.
%
%   See also ACLUST and ACLUSTRUN.

% Directories
global DATAPATH
inpdir1 = [DATAPATH 'LentiFenti\Cluster\mat\'];   % mat files
inpdir2 = [DATAPATH 'LentiFenti\Cluster\fig\'];   % fig files
resdir1 = [DATAPATH 'LentiFenti\Cluster\dec\'];   % clustercut ('dec')
resdir2 = [DATAPATH 'LentiFenti\Cluster\ac\'];    % aftercut ('ac')
resdir3 = [DATAPATH 'LentiFenti\Cluster\mat2\'];  % corrected burst param., based on aftercut
resdir4 = [DATAPATH 'LentiFenti\Cluster\xls\'];   % excel files
mm = pwd;

% Filelist
[files1 files_short1] = filelist(inpdir1);
[files2 files_short2] = filelist(inpdir2);
sf = length(files_short1);

% Cluster analysis
for o = 1:sf
    fname1 = files_short1{o};
    cd(inpdir1)
    load(fname1)     % load mat
    fname2 = files_short2{o};
    cd(inpdir2)
    open(fname2)     % load fig
    if ~isequal(fname1(1:end-4),fname2(1:end-4))
        error('Input mismatch.')
    else
        fname = fname1(1:end-10);
    end
    
    dec = [];
    while isempty(dec)
        dec = input('Number of clusters: ','s');
        dec = str2num(dec);
    end
    
    if ~isequal(dec,0)
        yn = 'y';
        while isequal(yn,'y')
            disp(['Burst number: ' num2str(length(After{dec}))])
            bno = [];
            while isempty(bno)
                bno = input('Number of bins: ','s');
                bno = str2num(bno);
            end
            figure
            hist(After{dec},bno)
            yn = [];
            while isempty(yn)
                yn = input('New histogram? [y/n] ','s');
                if ~(isequal(yn,'n') | isequal(yn,'y'))
                    yn = [];
                end
            end
        end

        ac = [];
        while isempty(ac)
            ac = input('Number of missed spikes: ','s');
            ac = str2num(ac);
        end
    else
        ac = 0;
    end
    
    if ~isequal(dec,0)
        burstparam(Burst{dec},vdisc,After{dec},ac,resdir3,resdir4,fname)    % burst mtx. correction
    end
    
% Save
    cd(resdir1)
    fn = [fname '_DEC.mat'];
    save(fn,'dec')
    
    cd(resdir2)
    fn = [fname '_AC.mat'];
    fn2 = [fname '_AC.fig'];
    save(fn,'ac')
    saveas(gcf,fn2)
    
    close all
end
cd(mm)

% -------------------------------------------------------------------------
function [files2, files2_short] = filelist(inpdir)

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
function burstparam(burst,vdisc,after,ac,resdir1,resdir2,fname)

% Corrigate burst matrix
isi = diff(vdisc);
afs = sort(after);
next = 1;
tburst = burst;
for k = 1:ac
    ci = round(afs(k)*20000);
    fni = find(isi==ci);
    if length(fni) > 1
        im = ismember(fni,tburst);
        fim = find(im);
        fni = fni(fim(next));
        if next == length(fim)
            next = 1;
%            tburst = burst;
        else
            next = next + 1;
        end
    else
%        tburst = burst;
    end
    vfn = vdisc(fni+1);
    vburst = vdisc(burst);
    fnb = find(vburst<vfn,1,'last');
    burst(fnb) = burst(fnb) + 1;
end
Burst = burst;

% Intraburst intervals
intraburstiv = [];
burstnum = size(burst,2);
intraburstnum = zeros(1,burstnum);
isimtx = {};
for j = 1:burstnum    %computing intraburstiv and allfirstspike
    b = vdisc(burst(1,j):burst(2,j));
    intraburstiv = [intraburstiv diff(b)];
    intraburstnum(j) = length(b);   %intraburst spike number
    for k = 1:length(b) - 1     % ISI matrix
        if size(isimtx,1) < (length(b) - 1) | size(isimtx,2) < k
            isimtx{length(b)-1,k} = [];
        end
        isimtx{length(b)-1,k} = [isimtx{length(b)-1,k} (vdisc(burst(1,j)+k)-vdisc(burst(1,j)+k-1))/20]; %in ms
    end
end
burstlength = (vdisc(burst(2,:)) - vdisc(burst(1,:))) / 20;  %in ms 
intraburstfreq = (intraburstnum - 1) ./ (burstlength / 1000);   %in Hz

% Burst frequency
if ~isequal(burstnum,0)
    firstspikefreq = 20000 * (burstnum - 1)  / (vdisc(burst(2,end)) - vdisc(burst(1,1)));   %in Hz
else
    firstspikefreq = NaN;
end

% Burst parameters
Burstiness = (length(intraburstiv) + burstnum) / length(vdisc);
IntraBurstFrequency.mean = mean(intraburstfreq);
IntraBurstFrequency.sd = std(intraburstfreq);
IntraBurstFrequency.all = intraburstfreq;
IntraBurstSpikeNumber.mean = mean(intraburstnum);
IntraBurstSpikeNumber.sd = std(intraburstnum);
IntraBurstSpikeNumber.all = intraburstnum;
BurstLength.mean = mean(burstlength);
BurstLength.sd = std(burstlength);
BurstLength.all = burstlength;
BurstFrequency = firstspikefreq;
for x = 1:size(isimtx,1)
    warning off
    for y = 1:size(isimtx,2)
        IsiMatrix.mean(x,y) = mean(isimtx{x,y});
        IsiMatrix.sd(x,y) = std(isimtx{x,y});
        IsiMatrix.num(x,y) = length(isimtx{x,y});
    end
    warning backtrace
end

% Save
mm = pwd;
cd(resdir1)
fn = [fname '_CLUST2.mat'];
save(fn,'Burstiness','IntraBurstFrequency','IntraBurstSpikeNumber','BurstLength',...
    'BurstFrequency','IsiMatrix','Burst','vdisc')

cd(resdir2)
cmps = strread(fname,'%s','delimiter','_');
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
cd(mm)