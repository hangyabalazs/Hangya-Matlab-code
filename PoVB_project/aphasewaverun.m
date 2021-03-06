function aphasewaverun
%APHASEWAVERUN    Runs APHASE, ASTANORM and WAVELET on a sequence of files.
%   APHASEWAVERUN loads EEG and discriminated unit, calls APHASE and ASTANORM 
%   and saves results. It also calculates and saves wavelet power spectrum.
%   Edit the program code to modify the input and result
%   directories.
%
%   See also APHASERUN, APHASE and ASTANORM.

% Directories
inpdir1 = 'X:\In_Vivo\_analysis\acsady_csoport\Ivan PO\mat\';
inpdir2 = 'd:\_analysis\matlab_data\Hajni\Disc\';
resdir1 = 'd:\_analysis\matlab_data\Hajni\Phase\';
resdir2 = 'd:\_analysis\matlab_data\Hajni\STA\';
resdir3 = 'd:\_analysis\matlab_data\Hajni\Wavelet\';
mm = pwd;

% Filelist
[files1 files_short1] = filelist(inpdir1);
[files2 files_short2] = filelist2(inpdir2);
files_short = intersect(files_short1,files_short2);
sf = length(files_short);

% Figures
H1 = figure;
H2 = figure;
H3 = figure;
H4 = figure;

% Main
for o = 1:sf
    fname = files_short{o};
    ff = [inpdir1 fname];
    load(ff)
    ff2 = [inpdir2 fname(1:end-4) '_d.mat'];
    load(ff2)
    eeg = data(:,2)';
    eeg = eeg(1:20:end);    % downsample on 1000 Hz
    clear data
    vdisc = round(vdisc/20);
    sr = 1000;
    
    [aang cyclen1 cyclen2 cl] = aphase(eeg,vdisc,sr);    % PHASE
    n = length(aang);
    ftm = sum(exp(1).^(i*aang)) / n;    % first trigonometric moment
    ang = angle(ftm);   % mean angle
    mvl = abs(ftm);     % mean resultant length
    aang = aang * 180 / pi;
    ang = ang * 180 / pi;
    [nm,xout] = hist(aang,18);   % phase histogram
    figure(H1);
    bar(xout,nm/length(aang))
    ach = allchild(H1);     % figure title
    ax = findobj(ach,'type','axes');
    cmps = strread(fname,'%s','delimiter','_');
    titlestr = [];
    for tt = 1:length(cmps)
        titlestr = [titlestr ' ' cmps{tt}];
    end
    title(ax(end),titlestr)
    x_lim = xlim;
    y_lim = ylim;
    axis([-200 200 y_lim(1) y_lim(2)])
    str = ['\it{Mean angle: }' '\bf ' num2str(ang)];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','red')
    str = ['\it{Mean resultant length: }' '\bf ' num2str(mvl)];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/12,str,'Color','red')
    str = ['\it{n: }' '\bf ' num2str(n)];
    text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/6,str,'Color','red')
    str = ['\it{Mean cycle length before: }' '\bf ' num2str(cyclen1)];
    text(-160,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','cyan')
    str = ['\it{Mean cycle length after: }' '\bf ' num2str(cyclen2)];
    text(-160,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/12,str,'Color','cyan')
    
    figure(H3)      % cycle discard statistics
    ccl = [cl{1} cl{2} cl{3} cl{4} cl{5} cl{6} cl{7}];
    lc1 = length(cl{1});
    lc2 = length(cl{2});
    lc3 = length(cl{3});
    lc4 = length(cl{4});
    lc5 = length(cl{5});
    lc6 = length(cl{6});
    lc7 = length(cl{7});
    gr = [zeros(1,lc1) ones(1,lc2) 2*ones(1,lc3) 3*ones(1,lc4)...
        4*ones(1,lc5) 5*ones(1,lc6) 6*ones(1,lc7)];
    str1 = ['Original(' num2str(lc1) ')'];
    str2 = ['Low EEG(' num2str(lc2) ')'];
    str3 = ['Non-monotonity(' num2str(lc3) ')'];
    str4 = ['Only non-mono.(' num2str(lc4) ')'];
    str5 = ['Short(' num2str(lc5) ')'];
    str6 = ['All discarded(' num2str(lc6) ')'];
    str7 = ['All remaining(' num2str(lc7) ')'];
    boxplot(ccl,gr,'labels',[{str1} {str2} {str3} {str4} {str5} {str6} {str7}])
    ach = allchild(H3);     % figure title
    ax = findobj(ach,'type','axes');
    title(ax(end),titlestr)
    
    wn = 2 * sr;    % 2 sec. window
    [sta sta_index1 sta_index2 nn] = astanorm(vdisc,eeg,wn);    % STA
    time = linspace(-wn/sr/2,wn/sr/2,length(sta));
    figure(H2);
    plot(time,sta)
    ach = allchild(H2);     % figure title
    ax = findobj(ach,'type','axes');
    title(ax(end),titlestr)
    x_lim = xlim;
    y_lim = ylim;
    str = ['\it{Max-mean: }' '\bf ' num2str(sta_index1)];
    text(0.1,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','red')
    str = ['\it{Max: }' '\bf ' num2str(sta_index2)];
    text(0.1,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/12,str,'Color','red')
    str = ['\it{n: }' '\bf ' num2str(nn)];
    text(0.1,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/6,str,'Color','red')
    
    [pow,phase,f] = eeg_wavelet(eeg);        % WAVELET
    figure(H4)
    imagesc(pow)
    b_rescaleaxis('Y',f)
    ach = allchild(H4);     % figure title
    ax = findobj(ach,'type','axes');
    title(ax(end),titlestr)
    
% Save
    cd(resdir1)
    fn11 = [fname(1:end-4) '_PHASE.jpg'];
    fn115 = [fname(1:end-4) '_PHASECRIT.jpg'];
    fn12 = [fname(1:end-4) '_PHASE.mat'];
    saveas(H1,fn11)
    saveas(H3,fn115)
    save(fn12,'aang','ang','mvl','cyclen1','cyclen2')
    
    cd(resdir2)
    fn21 = [fname(1:end-4) '_STA.jpg'];
    fn22 = [fname(1:end-4) '_STA.mat'];
    saveas(H2,fn21)
    save(fn22,'sta_index1','sta_index2')
    
    cd(resdir3)
    fn1 = [fname(1:end-4) '_WAVE.jpg'];
    saveas(H4,fn1)
end
close(H1)
close(H2)
close(H3)
close(H4)
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
function [pow,phase,f] = eeg_wavelet(dat)

% Prepare for wavelet transformation
variance = std(dat) ^ 2;
dat = (dat - mean(dat)) / sqrt(variance) ;
n = length(dat);
dt = 1 / 1000;
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
mif = 0.2;          %minimal intresting frequency
mis = find(f>mif);
mis = mis(end);     %maximal intristing scale
mother = 'Morlet';

% Wavelet transformation
[wave,period,scale,coi] = b_wavelet_new3(dat,dt,pad,dj,s0,j1,mother,param,mis);
pow = abs(wave).^2;
phase = angle(wave);