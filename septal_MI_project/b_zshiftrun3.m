function b_zshiftrun3
%ZSHIFTRUN3    Runs ZSHIFT on a set of files, calculates phase angles.
%   ZSHIFTRUN3 uses an input and a result directory: user can set them
%   through editing the program code.
%
%   Maximum localizations of Z arrays for theta and non-theta segments,
%   there means and standard deviations for one cell and along all cells
%   are getting saved in mat files. The plots are getting saved in jpg
%   files.
%   Hilbert and croosswavelet phase angles are also saved as well as 
%   correlation between zshift values and angles.
%
%   ZSHIFTRUN3 uses its own ZSHIFT instead of ZSHIFT2.
%
%   See also ZSHIFT2, ZSHIFTRUN and ZSHIFTRUN2.

% Input argument check
error(nargchk(0,0,nargin))

% Directories
global DATAPATH
global DATADIR
global DATADIR2
% where = DATADIR2;    %Here are the discriminated data files
where = 'F:\raw_data\disctemp\';
files = dir(where);
files = files(3:end);
sf = length(files);
mmm = pwd;
cd([DATAPATH,'Zshift5/']);  %Here are the results
create_subdir;         %create subdirectories

% Progress indicator
wb = waitbar(0,'Running ZSHIFTRUN3...','Position',[360 250 315 50]);    %Progress indicator
global WB
WB(end+1) = wb;

% Open figure
H = figure;

% Open text file
npt = input('Discard existing content or append data while writing text files? /discard:  ENTER, append: a/','s');
if isempty(npt)
    fid_th = fopen('zshift_theta.txt','w');
    fid_no = fopen('zshift_nontheta.txt','w');
elseif npt == 'a',
    fid_th = fopen('zshift_theta.txt','a');
    fid_no = fopen('zshift_nontheta.txt','a');
else
    error('Unexpected answer for input.')
end

% Import
list = {};      %initialize list of analysed files
ThetaMean = zeros(1,sf);
ThetaStd = zeros(1,sf);
CumThetaAng = [];
CumThetaZMaxLoc = [];
CumNonThetaAng = [];
CumNonThetaZMaxLoc = [];
for o = 1:sf    %"CELL CYCLE"
    fname = files(o).name;
    ffnm = [where fname];
    filenam = fname(1:6)
    pont = findstr(fname,'.');
    fn_noext = fname(1:pont(1)-1);
    load(ffnm);
    datinx1 = 1;      %first point of the interval
    datinx2 = length(eeg);       %last point of the interval
        
    islist = 1;     %initialize 'islist' for the decision of adding the cell name to the output list

% Loc. of first and last points of theta intervals
    fln = ['THETA_SEGMENTS_',filenam];
    ff = fullfile(DATAPATH,'Wavelet2\theta_segments',fln);
    load(ff)
    
    segments = ThetaSegments;
    if ~isempty(segments)
        segments = short_killer(segments);      % leave 200 - 200 ms edges then skip segments
    end

% Run 'zshift' on theta intervals
    th_index = size(segments,2);
    ThetaZMaxLoc = zeros(1,th_index);
    ThetaZMax = zeros(1,th_index);
    ThetaAng = zeros(2,th_index);
    ThetaMvl = zeros(2,th_index);
    theta_skiplist = {};
    for t = 1:th_index      % "THETA SEGMENT CYCLE"
        segfirst = segments(1,t);
        seglast = segments(2,t);
        eeg_input = eeg(segfirst:seglast);
        vdisc_input = vdisc(find(vdisc>=segfirst&vdisc<=seglast)) - segfirst + 1;
        
        if length(vdisc_input) ~= 0
            [zmaxloc,zmax,hang,wang,hmvl,wmvl] = zshift(H,eeg_input,vdisc_input,filenam,segfirst,seglast);
        end
        ThetaZMaxLoc(t) = zmaxloc;
        ThetaZMax(t) = zmax;
        CumThetaZMaxLoc = [CumThetaZMaxLoc zmaxloc];
        ThetaAng(1,t) = hang;
        ThetaAng(2,t) = wang;
        CumThetaAng = [CumThetaAng hang];
        ThetaMvl(1,t) = hmvl;
        ThetaMvl(2,t) = wmvl;
        cd jpg
        eval(['saveas(H,''THETA_',filenam,'_',num2str(segfirst),'_',num2str(seglast),'_ZSHIFT.jpg'')']);
        cd ..
        if zmaxloc == NaN   % add skipedd file to a list
            theta_skiplist{end+1} = [filenam,'_',num2str(segfirst),'_',num2str(seglast)];
        end
        fprintf(fid_th,'%s %f %f %f %f %f %f %f %f \n',...
        filenam,segfirst,seglast,...
        zmaxloc,zmax,hang,hmvl,wang,wmvl);      % write text file
    end     % end of "theta segment cycle"
    if ~isempty(ThetaZMaxLoc(find(~isnan(ThetaZMaxLoc))))
        ThetaMean(o) = b_mean_nonnan(ThetaZMaxLoc);
        ThetaStd(o) = b_std_nonnan(ThetaZMaxLoc);
    else
        ThetaMean(o) = NaN;
        ThetaStd(o) = NaN;
    end
    
% Loc. of first and last points of non-theta intervals
    fln = ['NONTHETA_SEGMENTS_',filenam];
    ff = fullfile(DATAPATH,'Wavelet2\nontheta_segments',fln);
    load(ff)
    
    segments = NonThetaSegments;
    if ~isempty(segments)
        segments = short_killer(segments);      % leave 200 - 200 ms edges then skip segments
    end

% Run 'zshift' on non-theta intervals
    no_index = size(segments,2);
    NonThetaZMaxLoc = zeros(1,no_index);
    NonThetaZMax = zeros(1,no_index);
    NonThetaAng = zeros(2,no_index);
    NonThetaMvl = zeros(2,no_index);
    nontheta_skiplist = {};
    for t = 1:no_index      % "THETA SEGMENT CYCLE"
        segfirst = segments(1,t);
        seglast = segments(2,t);
        eeg_input = eeg(segfirst:seglast);
        vdisc_input = vdisc(find(vdisc>=segfirst&vdisc<=seglast)) - segfirst + 1;
        
        if length(vdisc_input) ~= 0
            [zmaxloc,zmax,hang,wang,hmvl,wmvl] = zshift(H,eeg_input,vdisc_input,filenam,segfirst,seglast);
        end
        NonThetaZMaxLoc(t) = zmaxloc;
        NonThetaZMax(t) = zmax;
        CumNonThetaZMaxLoc = [CumNonThetaZMaxLoc zmaxloc];
        NonThetaAng(1,t) = hang;
        NonThetaAng(2,t) = wang;
        CumNonThetaAng = [CumNonThetaAng hang];
        NonThetaMvl(1,t) = hmvl;
        NonThetaMvl(2,t) = wmvl;
        cd jpg
        eval(['saveas(H,''NONTHETA_',filenam,'_',num2str(segfirst),'_',num2str(seglast),'_ZSHIFT.jpg'')']);
        cd ..
        if zmaxloc == NaN   % add skipedd file to a list
            nontheta_skiplist{end+1} = [filenam,'_',num2str(segfirst),'_',num2str(seglast)];
        end
        fprintf(fid_no,'%s %f %f %f %f %f %f %f %f \n',...
        filenam,segfirst,seglast,...
        zmaxloc,zmax,hang,hmvl,wang,wmvl);      % write text file
    end     % end of "theta segment cycle"
    if ~isempty(NonThetaZMaxLoc(find(~isnan(NonThetaZMaxLoc))))
        NonThetaMean(o) = b_mean_nonnan(NonThetaZMaxLoc);
        NonThetaStd(o) = b_std_nonnan(NonThetaZMaxLoc);
    else
        NonThetaMean(o) = NaN;
        NonThetaStd(o) = NaN;
    end
    
    cd mat
    str = [filenam '_ZSHIFT'];
    save(str,'ThetaZMaxLoc','NonThetaZMaxLoc','ThetaZMax','NonThetaZMax',...
        'ThetaAng','NonThetaAng','ThetaMvl','NonThetaMvl')
    str = [filenam '_ANG'];
    save(str,'ThetaAng','NonThetaAng','ThetaMvl','NonThetaMvl')
    cd ..
    
    if islist    %generate list of analysed files
        list{end+1} = filenam;
    end
    waitbar(o/sf)
end

AllThetaMean = b_mean_nonnan(ThetaMean);
AllThetaStd = b_std_nonnan(ThetaMean);
AllNonThetaMean = b_mean_nonnan(NonThetaMean);
AllNonThetaStd = b_std_nonnan(NonThetaMean);
cd mat
str = 'all_ZSHIFT';
save(str,'ThetaMean','ThetaStd','NonThetaMean','NonThetaStd','AllThetaMean','AllThetaStd',...
    'AllNonThetaMean','AllNonThetaStd')
save ThetaSkipList theta_skiplist;
save NonThetaSkipList nontheta_skiplist;

% Correlation between phase angles and 'zmaxloc'
fn1 = find(isnan(CumThetaAng));
CumThetaAng(fn1) = [];
CumThetaZMaxLoc(fn1) = [];
fn2 = find(isnan(CumThetaZMaxLoc));
CumThetaAng(fn2) = [];
CumThetaZMaxLoc(fn2) = [];
ThetaAngZCorr = corr(CumThetaAng',CumThetaZMaxLoc');

fn1 = find(isnan(CumNonThetaAng));
CumNonThetaAng(fn1) = [];
CumNonThetaZMaxLoc(fn1) = [];
fn2 = find(isnan(CumNonThetaZMaxLoc));
CumNonThetaAng(fn2) = [];
CumNonThetaZMaxLoc(fn2) = [];
NonThetaAngZCorr = corr(CumNonThetaAng',CumNonThetaZMaxLoc');
save CorrData ThetaAngZCorr NonThetaAngZCorr CumThetaAng CumNonThetaAng...
    CumThetaZMaxLoc CumNonThetaZMaxLoc
cd ..

% Finalize
datum = datevec(datestr(now));
save list list datum
cd(mmm)
fclose(fid_th);
fclose(fid_no);
close(wb)
close(H)



% ----------------------------------------------------------------------------------
function create_subdir

% Create subdirectories
if ~b_isdir2('jpg')
    mkdir jpg
end
if ~b_isdir2('mat')
    mkdir mat
end



% ----------------------------------------------------------------------------------
function segments = short_killer(segments)

% Skip short segments
int = segments;
int1 = int(1,:) + 200;    % leaving 200 - 200 ms edges of each segment
int2 = int(2,:) - 200;
difint = int2 - int1;
fd = find(difint<50000);         % leaving segments shorter than 5 sec.
int1(fd) = [];
int2(fd) = [];
segments = [int1; int2];



% ----------------------------------------------------------------------------------
function [zmaxloc,zmax,hang,wang,hmvl,wmvl] = zshift(H,eeg,vdisc,filenam,segfirst,seglast)

% Filtering EEG
flt = fir1(512,[2/5000 8/5000]);
feeg = filtfilt(flt,1,eeg);

% Hilbert transformation of the eeg
ahee = angle(hilbert(feeg));
aheedeg = ahee * (180 / pi);

% Phase angles - Hilbert
bahee = ahee(vdisc);
n = length(bahee);
ftm = sum(exp(1).^(i*bahee)) / n;    % first trigonometric moment
hang = angle(ftm);   % mean angle
hmvl = abs(ftm);     % mean resultant length

% Phase angles - Crosswavelet
[wave_eeg,f] = eeg_wavelet(eeg(1:50:end));
[wave_unit,f] = unit_wavelet(vdisc,length(eeg));
if size(wave_eeg,2) < size(wave_unit,2)
    wave_unit = wave_unit(:,1:end-1);
end
if size(wave_unit,2) < size(wave_eeg,2)
    wave_eeg = wave_eeg(:,1:end-1);
end
wave_cross = wave_eeg .* conj(wave_unit);

fnd = find(f>6);    % find theta band
pwind1 = fnd(end);
fnd = find(f<2.5);
pwind2 = fnd(1);

wph = angle(wave_cross(pwind1:pwind2,:));
wphh = reshape(wph,1,numel(wph));

n = length(wphh);
ftm = sum(exp(1).^(i*wphh)) / n;    % first trigonometric moment
wang = angle(ftm);   % mean angle
wmvl = abs(ftm);     % mean resultant length

% Shift unit
Z = [];
p = [];
T = [-10000:10:10000];
for t = T
    vt = vdisc(find(vdisc+t>0&vdisc+t<length(eeg))) + t;  % shift unit
    if isempty(vt)      % Skip segment, if all spikes fall out of range after unit shift
        zmaxloc = NaN;
        zmax = NaN;
        figure(H);
        plot(0,0)
        text(0,0,'All spikes out of range after unit shift.','Color','red',...
            'HorizontalAlignment','center');
        return
    end
    bang = ahee(vt);

% Mean resultant length
    n = length(bang);
    ftm = sum(exp(1).^(i*bang)) / n;    % first trigonometric moment
    mrl = abs(ftm);     % mean resultant length
    z = n * (mrl ^ 2);  % Rayleigh's Z statistic
    Z(end+1) = z;
    p(end+1) = exp(1) ^ (-1 * z) * (1 + (2 * z - z ^ 2) / ...
        (4 * n) - (24 * z - 132 * z ^ 2 + 76 * z ^ 3 - 9 * z ^ 4) / (288 * n ^ 2));
end

% Plot result
figure(H);
plot(T,Z)

sl = 0.05;  % level of significance
siglev = -1 * log(sl);
hold on
line(xlim,[siglev siglev],'Color','k','LineStyle','--');

zmax = max(Z);
zmxlc = find(Z==zmax);
x_lim = xlim;
y_lim = ylim;
text(x_lim(1)+2*(x_lim(2)-x_lim(1))/3,y_lim(1)+2*(y_lim(2)-y_lim(1))/3,num2str(T(zmxlc(1))));
tt = [filenam(1:3) ' ' filenam(5:6) ' ' num2str(segfirst) ' ' num2str(seglast)];
title(tt)
hold off

% Skip file, if phase connection is not significant at the maximum localization
if Z(zmxlc) < siglev
    zmaxloc = NaN;
else
    zmaxloc = T(zmxlc);
end



% --------------------------------------------------------------------------------------------
function [wave,f] = eeg_wavelet(dat)

% Prepare for wavelet transformation
variance = std(dat) ^ 2;
dat = (dat - mean(dat)) / sqrt(variance) ;
n = length(dat);
dt = 1 / 200;
pad = 1;
dj = 0.1;    
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



% --------------------------------------------------------------------------------------------
function [wave,f] = unit_wavelet(vdisc,lenu)

% Sinc convolution
fs = 10000;     % unit
dto = 1 / fs;
fcut = 100; 
fsnew = 200;
dtnew = 1 / fsnew;
fsold = 10000;
fsratio = fsnew / fsold;
told = vdisc * dto * fcut;
tnew = (1:lenu*fsratio) * dtnew * fcut;
lentold = length(told);
zint = 0;
for i = 1:lentold
    zint = zint + sinc(tnew-told(i));
end

% Prepare for wavelet transformation
n = length(zint);
dt = 1 / fsnew;
pad = 1;
dj = 0.1;
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
[wave,period,scale,coi] = b_wavelet_new3(zint,dt,pad,dj,s0,j1,mother,param,mis);