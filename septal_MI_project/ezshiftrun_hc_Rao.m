function ezshiftrun_hc_Rao
%EZSHIFTRUN_HC_RAO    Runs ZSHIFT on a set of files, calculates phase angles.
%   EZSHIFTRUN_HC_RAO uses an input and a result directory: user can set them
%   through editing the program code.
%
%   Maximum localizations of Z arrays for theta and non-theta segments,
%   their means and standard deviations for one cell and along all cells
%   are getting saved in mat files. The plots are getting saved in jpg
%   files.
%   Hilbert phase angles are also saved as well as regression statistics
%   between zshift values and angles/estimated concentration parameters.
%
%   EZSHIFTRUN_HC_RAO uses its own ZSHIFT instead of ZSHIFT2.
%
%   EZSHIFTRUN_HC_RAO uses Rao's spacing test (p = 0.001).
%
%   See also ZSHIFT2 and EZSHIFTRUN.

% Input argument check
error(nargchk(0,0,nargin))
dbstop if error

% Directories
global DATAPATH
global DATADIR
global DATADIR2
where = [DATADIR 'mi_zshift_data\discriminated_hc\'];    %Here are the discriminated data files
files = dir(where);
files = files(3:end);
sf = length(files);
mmm = pwd;
cd([DATAPATH,'Ezshift_hc_Rao2\']);  %Here are the results
create_subdir;         %create subdirectories

% Progress indicator
wb = waitbar(0,'Running EZSHIFTRUN HC RAO...','Position',[360 250 315 50]);    %Progress indicator
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
xlsout_no = cell(1,7);     % initialize Excel output
xlsout_no(1,1:7) = [{'Filename'} {'first point'} {'last point'} {'Z-shift'}...
    {'Z max'} {'angle'} {'mvl'}];
xlsout_th = cell(1,9);
xlsout_th(1,1:7) = [{'Filename'} {'first point'} {'last point'} {'Z-shift'}...
    {'Z max'} {'angle'} {'mvl'}];

% Import
list = {};      %initialize list of analysed files
ThetaMean = zeros(1,sf);
ThetaStd = zeros(1,sf);
CumThetaAng = [];
CumThetaZMaxLoc = [];
CumThetaKappa = [];
CumNonThetaAng = [];
CumNonThetaZMaxLoc = [];
for o = 1:sf    %"CELL CYCLE"
    fname = files(o).name;
    ffnm = [where fname];
    filenam = fname(1:6)
    filenam(findstr(filenam,'H')) = '_';   % 3ch cellnames contain 'H'
    pont = findstr(fname,'.');
    fn_noext = fname(1:pont(1)-1);
    load(ffnm);
    datinx1 = 1;      %first point of the interval
    datinx2 = length(eeg);       %last point of the interval
        
    islist = 1;     %initialize 'islist' for the decision of adding the cell name to the output list
    
% Loc. of first and last points of theta intervals
    fln = ['THETA_SEGMENTS_',filenam];
    ff = fullfile(DATAPATH,'Wavelet_hc\theta_segments',fln);
    load(ff)
    
    segments = ThetaSegments;
    if ~isempty(segments)
        segments = short_killer(segments);      % leave 200 - 200 ms edges then skip segments
    end

% Run 'zshift' on theta intervals
    th_index = size(segments,2);
    ThetaZMaxLoc = zeros(1,th_index);
    ThetaZMax = zeros(1,th_index);
    ThetaAng = zeros(1,th_index);
    ThetaMvl = zeros(1,th_index);
    ThetaKappa = zeros(1,th_index);
    theta_skiplist = {};
    for t = 1:th_index      % "THETA SEGMENT CYCLE"
        segfirst = segments(1,t);
        seglast = segments(2,t);
        eeg_input = eeg(segfirst:seglast);
        vdisc_input = vdisc(find(vdisc>=segfirst-10000&vdisc<=seglast+10000)) - segfirst + 1;
        
        if length(vdisc_input) ~= 0
            [zmaxloc,zmax,hang,hmvl] = zshift(H,eeg_input,vdisc_input,filenam,segfirst,seglast);
        end
        ThetaZMaxLoc(t) = zmaxloc;
        ThetaZMax(t) = zmax;
        CumThetaZMaxLoc = [CumThetaZMaxLoc zmaxloc];
        ThetaAng(1,t) = hang;
        CumThetaAng = [CumThetaAng hang];
        ThetaMvl(1,t) = hmvl;
        if hmvl < 0.53
            ThetaKappa(1,t) = 2 * hmvl + hmvl^3 + 5 * hmvl^5 / 6;
        elseif hmvl < 0.85
            ThetaKappa(1,t) = -0.4 + 1.39 * hmvl + 0.43 / (1 - hmvl);
        else
            ThetaKappa(1,t) = 1 / (hmvl^3 - 4 * hmvl^2 + 3 * hmvl);
        end
        CumThetaKappa = [CumThetaKappa ThetaKappa(1,t)];
        cd jpg
        dbclear if error
        eval(['saveas(H,''THETA_',filenam,'_',num2str(segfirst),'_',num2str(seglast),'_ZSHIFT.jpg'')']);
        dbstop if error
        cd ..
        if zmaxloc == NaN   % add skipedd file to a list
            theta_skiplist{end+1} = [filenam,'_',num2str(segfirst),'_',num2str(seglast)];
        end
        fprintf(fid_th,'%s %f %f %f %f %f %f \n',...
            filenam,segfirst,seglast,...
            zmaxloc,zmax,hang,hmvl);      % write text file
        xlsout_th(end+1,1:7) = [{filenam} {segfirst} {seglast} {zmaxloc}...
            {zmax} {hang} {hmvl}];
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
    ff = fullfile(DATAPATH,'Wavelet_hc\nontheta_segments',fln);
    load(ff)
    
    segments = NonThetaSegments;
    if ~isempty(segments)
        segments = short_killer(segments);      % leave 200 - 200 ms edges then skip segments
    end

% Run 'zshift' on non-theta intervals
    no_index = size(segments,2);
    NonThetaZMaxLoc = zeros(1,no_index);
    NonThetaZMax = zeros(1,no_index);
    NonThetaAng = zeros(1,no_index);
    NonThetaMvl = zeros(1,no_index);
    nontheta_skiplist = {};
    for t = 1:no_index      % "THETA SEGMENT CYCLE"
        segfirst = segments(1,t);
        seglast = segments(2,t);
        eeg_input = eeg(segfirst:seglast);
        vdisc_input = vdisc(find(vdisc>=segfirst-10000&vdisc<=seglast+10000)) - segfirst + 1;
        
        if length(vdisc_input) ~= 0
            [zmaxloc,zmax,hang,hmvl] = zshift(H,eeg_input,vdisc_input,filenam,segfirst,seglast);
        end
        NonThetaZMaxLoc(t) = zmaxloc;
        NonThetaZMax(t) = zmax;
        CumNonThetaZMaxLoc = [CumNonThetaZMaxLoc zmaxloc];
        NonThetaAng(1,t) = hang;
        CumNonThetaAng = [CumNonThetaAng hang];
        NonThetaMvl(1,t) = hmvl;
        cd jpg
        dbclear if error
        eval(['saveas(H,''NONTHETA_',filenam,'_',num2str(segfirst),'_',num2str(seglast),'_ZSHIFT.jpg'')']);
        dbstop if error
        cd ..
        if zmaxloc == NaN   % add skipedd file to a list
            nontheta_skiplist{end+1} = [filenam,'_',num2str(segfirst),'_',num2str(seglast)];
        end
        fprintf(fid_no,'%s %f %f %f %f %f %f \n',...
            filenam,segfirst,seglast,...
            zmaxloc,zmax,hang,hmvl);      % write text file
        xlsout_no(end+1,1:7) = [{filenam} {segfirst} {seglast} {zmaxloc}...
            {zmax} {hang} {hmvl}];
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
CumThetaKappa(fn1) = [];
fn2 = find(isnan(CumThetaZMaxLoc));
CumThetaAng(fn2) = [];
CumThetaZMaxLoc(fn2) = [];
CumThetaKappa(fn2) = [];
[b,bint,r,rint,stats] = regress(CumThetaZMaxLoc',[ones(length(CumThetaZMaxLoc),1) sin(CumThetaAng') cos(CumThetaAng')]);
ThetaAngZReg = stats;
[b,bint,r,rint,stats] = regress(CumThetaZMaxLoc',[ones(length(CumThetaZMaxLoc),1),CumThetaKappa']);
ThetaKappaZReg = stats;
HS1 = figure;
plot(CumThetaZMaxLoc,CumThetaAng,'.')
dbclear if error
saveas(HS1,'zshift_vs_ang')
HS2 = figure;
plot(CumThetaZMaxLoc,CumThetaKappa,'.')
saveas(HS2,'zshift_vs_kappa')
dbstop if error

fn1 = find(isnan(CumNonThetaAng));
CumNonThetaAng(fn1) = [];
CumNonThetaZMaxLoc(fn1) = [];
fn2 = find(isnan(CumNonThetaZMaxLoc));
CumNonThetaAng(fn2) = [];
CumNonThetaZMaxLoc(fn2) = [];
[b,bint,r,rint,stats] = regress(CumNonThetaZMaxLoc',[ones(length(CumNonThetaZMaxLoc),1) sin(CumNonThetaAng') cos(CumNonThetaAng')]);
NonThetaAngZReg = stats;
save RegData ThetaAngZReg ThetaKappaZReg NonThetaAngZReg CumThetaAng CumNonThetaAng...
    CumThetaZMaxLoc CumNonThetaZMaxLoc
xlswrite('zshift_theta.xls',xlsout_th)
xlswrite('zshift_noth.xls',xlsout_no)
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
function [zmaxloc,zmax,hang,hmvl] = zshift(H,eeg,vdisc,filenam,segfirst,seglast)

% Filtering EEG
flt = fir1(512,[2/5000 8/5000]);
feeg = filtfilt(flt,1,(eeg-mean(eeg))/std(eeg));

% Hilbert transformation of the eeg
ahee = angle(hilbert(feeg));
aheedeg = ahee * (180 / pi);

% Phase angles - Hilbert
bahee = ahee(vdisc(find(vdisc>0&vdisc<length(eeg))));
n = length(bahee);
ftm = sum(exp(1).^(i*bahee)) / n;    % first trigonometric moment
hang = angle(ftm);   % mean angle
hmvl = abs(ftm);     % mean resultant length

% Shift unit
U = [];
p = [];
T = [-10000:10:10000];
for t = T
    vt = vdisc(find(vdisc+t>0&vdisc+t<length(eeg))) + t;  % shift unit
    if isempty(vt)      % Skip segment, if all spikes fall out of range after unit shift
        zmaxloc = NaN;
        zmax = NaN;
        plot(0,0)
        text(0,0,'All spikes out of range after unit shift.','Color','red',...
            'HorizontalAlignment','center');
        return
    end
    bang = ahee(vt);

% Rao spacing test
    [U1,p1,p001] = rao(bang);
    U(end+1) = U1;
    p(end+1,1:2) = p1;
end

% Plot result
plot(T,U)

siglev = p001;  % level of significance: p=0.05
hold on
line(xlim,[siglev siglev],'Color','k','LineStyle','--');

zmax = max(U);
zmxlc = find(U==zmax);
x_lim = xlim;
y_lim = ylim;
try
    text(x_lim(1)+2*(x_lim(2)-x_lim(1))/3,y_lim(1)+2*(y_lim(2)-y_lim(1))/3,num2str(T(zmxlc(1))));
end
tt = [filenam(1:3) ' ' filenam(5:6) ' ' num2str(segfirst) ' ' num2str(seglast)];
title(tt)
hold off

% Skip file, if phase connection is not significant at the maximum localization
if isempty(zmxlc)
    zmaxloc = NaN;
elseif U(zmxlc) < siglev
    zmaxloc = NaN;
else
    zmaxloc = T(zmxlc);
end



% -------------------------------------------------------------------------
function [U,p,p001] = rao(x)
%RAO    Rao's Spacing Test.
%   [U,P,P001] = RAO(X) calculates Rao's U-Statistic (U) with the corresponding
%   p-value (P). P001 gives U at p=0.001 significance level. Input argument X
%   is the vector of phase angles in radians.

% Input argument check
error(nargchk(1,1,nargin))

% Angle
bang = mod(x,2*pi);     % from 0 to 2*pi (from -pi to pi: angle(cos(x)+i*sin(x)))

% Rao's Spacing Test
xdeg = x * 180 / pi;
xdegs = sort(mod(xdeg,360));
n = length(xdegs);
spacings = [diff(xdegs), xdegs(1)-xdegs(n)+360];
U = 0.5 * sum(abs(spacings-360/n));     % Test Statistic
if n < 4
    warning('Too small sample size.')
    U = NaN;
    p = NaN;
    p001 = NaN;
    return
elseif n <= 30
    trow = n - 3;
elseif n <= 32
    trow = 27;
elseif n <= 37
	trow = 28;
elseif n <= 42
	trow = 29;
elseif n <= 47
	trow = 30;
elseif n <= 62
	trow = 31;
elseif n <= 87
	trow = 32;
elseif n <= 125
	trow = 33;
elseif n <= 175
	trow = 34;
elseif n <= 250
	trow = 35;
elseif n <= 350
	trow = 36;
elseif n <= 450
	trow = 37;
elseif n <= 550
	trow = 38;
elseif n <= 650
	trow = 39;
elseif n <= 750
	trow = 40;
elseif n <= 850
	trow = 41;
elseif n <= 950
	trow = 42;
else
    trow = 43;
end
rao_table = raotable;
if U > rao_table(trow,1)
	p(1) = 0;
    p(2) = 0.001;
elseif U > rao_table(trow,2)
	p(1) = 0.001;
    p(2) = 0.01;
elseif U > rao_table(trow,3)
	p(1) = 0.01;
    p(2) = 0.05;
elseif U > rao_table(trow,4)
	p(1) = 0.05;
    p(2) = 0.10;
else
    p(1) = 0.10;
    p(2) = 1;
end
p001 = rao_table(trow,2);

% -------------------------------------------------------------------------
function RT = raotable
RT = [247.32, 221.14, 186.45, 168.02;...
    245.19, 211.93, 183.44, 168.66;...
    236.81, 206.79, 180.65, 166.30;...
    229.46, 202.55, 177.83, 165.05;...
    224.41, 198.46, 175.68, 163.56;...
    219.52, 195.27, 173.68, 162.36;...
    215.44, 192.37, 171.98, 161.23;...
    211.87, 189.88, 170.45, 160.24;...
    208.69, 187.66, 169.09, 159.33;...
    205.87, 185.68, 167.87, 158.50;...
    203.33, 183.90, 166.76, 157.75;...
    201.04, 182.28, 165.75, 157.06;...
    198.96, 180.81, 164.83, 156.43;...
    197.05, 179.46, 163.98, 155.84;...
    195.29, 178.22, 163.20, 155.29;...
    193.67, 177.08, 162.47, 154.78;...
    192.17, 176.01, 161.79, 154.31;...
    190.78, 175.02, 161.16, 153.86;...
    189.47, 174.10, 160.56, 153.44;...
    188.25, 173.23, 160.01, 153.05;...
    187.11, 172.41, 159.48, 152.68;...
    186.03, 171.64, 158.99, 152.32;...
    185.01, 170.92, 158.52, 151.99;...
    184.05, 170.23, 158.07, 151.67;...
    183.14, 169.58, 157.65, 151.37;...
    182.28, 168.96, 157.25, 151.08;...
    181.45, 168.38, 156.87, 150.80;...
    177.88, 165.81, 155.19, 149.59;...
    174.99, 163.73, 153.82, 148.60;...
    172.58, 162.00, 152.68, 147.76;...
    170.54, 160.53, 151.70, 147.05;...
    163.60, 155.49, 148.34, 144.56;...
    159.45, 152.46, 146.29, 143.03;...
    154.51, 148.84, 143.83, 141.18;...
    151.56, 146.67, 142.35, 140.06;...
    148.06, 144.09, 140.57, 138.71;...
    145.96, 142.54, 139.50, 137.89;...
    144.54, 141.48, 138.77, 137.33;...
    143.48, 140.70, 138.23, 136.91;...
    142.66, 140.09, 137.80, 136.59;...
    142.00, 139.60, 137.46, 136.33;...
    141.45, 139.19, 137.18, 136.11;...
    140.99, 138.84, 136.94, 135.92];