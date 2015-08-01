function axcorrisihrun
%AXCORRISIHRUN    Runs XCORR and calculates ISI hist. on a sequence of files.
%   AXCORRISIHRUN loads EEG and discriminated unit, calculates autocorrelation
%   and ISI histogram and saves results. Edit the program code to modify the
%   input and result directories.
%
%   See also XCORR.

% Directories
inpdir1 = 'X:\In_Vivo\_analysis\acsady_csoport\auj_istvan\mat_ket_xyl\';
inpdir2 = 'd:\_analysis\matlab_data\Andi\Disc\ketxyl\';
resdir2 = 'd:\_analysis\matlab_data\Andi\Ketxyl\ISIh\';
resdir1 = 'd:\_analysis\matlab_data\Andi\Ketxyl\Acorr\';
mm = pwd;

% Filelist
[files1 files_short1] = filelist(inpdir1);
[files2 files_short2] = filelist2(inpdir2);
files_short = intersect(files_short1,files_short2);
sf = length(files_short);

% Figures
H1 = figure;
H2 = figure;

% Main
for o = 1:sf
    fname = files_short{o};
    ff = [inpdir1 fname];
    load(ff)
    ff2 = [inpdir2 fname(1:end-4) '_d.mat'];
    load(ff2)
    lend = ceil(length(data)/10);     % downsample on 2000 Hz
    clear data
    vdisc = round(vdisc/10);
    sr = 2000;
    
    zunit = zeros(1,lend);
    zunit(vdisc) = 1;
    acr = xcorr(zunit,1*sr);
    acr(length(acr)/2+0.5) = [];
    acr = reshape(acr,length(acr)/200,200);
    sacr = sum(acr);
    figure(H1);
    bar(sacr)
    
    ach = allchild(H1);     % figure title
    ax = findobj(ach,'type','axes');
    cmps = strread(fname,'%s','delimiter','_');
    titlestr = [];
    for tt = 1:length(cmps)
        titlestr = [titlestr ' ' cmps{tt}];
    end
%     title(ax(end),titlestr)
    axis off
    
    isi = diff(vdisc);
    bno = fix(exp(0.626+0.4*log(length(isi)-1)));     % number of bins
    [n,isih] = hist(isi,bno);
    figure(H2);
    bar(isih,n)
    ach = allchild(H2);     % figure title
    ax = findobj(ach,'type','axes');
%     title(ax(end),titlestr)
    axis off
    
% Save
    cd(resdir1)
    fn11 = [fname(1:end-4) '_ACORR.fig'];
    fn12 = [fname(1:end-4) '_ACORR.jpg'];
    saveas(H1,fn11)
    saveas(H1,fn12)
    
    cd(resdir2)
    fn21 = [fname(1:end-4) '_ISIH.fig'];
    fn22 = [fname(1:end-4) '_ISIH.jpg'];
    saveas(H2,fn21)
    saveas(H2,fn22)
end
close(H1)
close(H2)
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