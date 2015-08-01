function azshiftrun5
%ZSHIFTRUN5    Runs ZSHIFT on a set of files, calculates phase angles.
%   ZSHIFTRUN5 uses an input and a result directory: user can set them
%   through editing the program code.
%
%   Maximum localizations of Z arrays for theta and non-theta segments,
%   there means and standard deviations for one cell and along all cells
%   are getting saved in mat files. The plots are getting saved in jpg
%   files.
%   Hilbert phase angles are also saved as well as regression statistics
%   between zshift values and angles/estimated concentration parameters.
%
%   ZSHIFTRUN5 uses its own ZSHIFT instead of ZSHIFT2.
%
%   ZSHIFTRUN5 uses Rayleigh's Z-test (p = 0.005).
%
%   See also ZSHIFT2, ZSHIFTRUN4 and ZSHIFTRUN_ONSET.

% Input argument check
error(nargchk(0,0,nargin))

% Directories
global DATAPATH
global DATADIR
where = [DATAPATH 'Andi\Disc\ketxyl\control\'];    %Here are the discriminated data files
% inpdir_eeg = 'X:\In_Vivo\_analysis\acsady_csoport\auj_istvan\mat_ket_xyl\';
inpdir_eeg = 'Y:\_Projects\AUJ_ISTVAN\DATA\MAT\control\';
files = b_filelist(where);
sf = length(files);
mmm = pwd;
cd([DATAPATH,'Andi\Ketxyl\Zshift\control\']);  %Here are the results
create_subdir;         %create subdirectories

% Progress indicator
wb = waitbar(0,'Running AZSHIFTRUN5...','Position',[360 250 315 50]);    %Progress indicator
global WB
WB(end+1) = wb;

% Open figure
H = figure;

% Open text file
npt = input('Discard existing content or append data while writing text files? /discard:  ENTER, append: a/','s');
if isempty(npt)
    fid = fopen('zshift.txt','w');
elseif npt == 'a',
    fid = fopen('zshift.txt','a');
else
    error('Unexpected answer for input.')
end

% Import
sr = 20000;
ZMean = zeros(1,sf);
ZStd = zeros(1,sf);
CumAng = [];
CumZMaxLoc = [];
CumKappa = [];
for o = 1:sf    %"CELL CYCLE"
    fname = files(o).name;
    ffnm = [where fname];
    cmps = strread(fname,'%s','delimiter','_');
    if length(cmps) == 1
        filenam = [cmps{1}];
        titlestr = [cmps{1}];
    elseif length(cmps) == 2
        filenam = [cmps{1} '_' cmps{2}];
        titlestr = [cmps{1} ' ' cmps{2}];
    else
        filenam = [cmps{1} '_' cmps{2} '_' cmps{3}];
        titlestr = [cmps{1} ' ' cmps{2} ' ' cmps{3}];
    end    
    pont = findstr(fname,'.');
    fn_noext = fname(1:pont(1)-1);
    load(ffnm);
    ff = [inpdir_eeg fname(1:end-6) '.mat'];       % load
    try
        load(ff)
    catch
        continue
    end
    eeg = data(:,2)';
    len = length(eeg);
    datinx1 = 1;      %first point of the interval
    datinx2 = length(eeg);       %last point of the interval
    
% Run 'zshift'
    seglen = 30 * sr;        % 30 sec. long segments
    lenr = floor(len/seglen);       % preallocation
    ZMaxLoc = zeros(1,lenr);    % preallocation
    ZMax = zeros(1,lenr);
    Ang = zeros(1,lenr);
    Mvl = zeros(1,lenr);
    Kappa = zeros(1,lenr);
    ind1 = [1:seglen:len];
    ind2 = ind1 + seglen -1;
    for t = 1:lenr      % "SEGMENT CYCLE"
        segfirst = ind1(t);
        seglast = ind2(t);
        eeg_input = eeg(segfirst:seglast);
        vdisc_input = vdisc(find(vdisc>=segfirst-1*sr&vdisc<=seglast+1*sr)) - segfirst + 1;
        
        if length(vdisc_input) ~= 0
            [zmaxloc,zmax,hang,hmvl] = zshift(H,eeg_input,vdisc_input,filenam,titlestr,segfirst,seglast,sr);
        end
        ZMaxLoc(t) = zmaxloc;
        ZMax(t) = zmax;
        CumZMaxLoc = [CumZMaxLoc zmaxloc];
        Ang(1,t) = hang;
        CumAng = [CumAng hang];
        Mvl(1,t) = hmvl;
        if hmvl < 0.53
            Kappa(1,t) = 2 * hmvl + hmvl^3 + 5 * hmvl^5 / 6;
        elseif hmvl < 0.85
            Kappa(1,t) = -0.4 + 1.39 * hmvl + 0.43 / (1 - hmvl);
        else
            Kappa(1,t) = 1 / (hmvl^3 - 4 * hmvl^2 + 3 * hmvl);
        end
        CumKappa = [CumKappa Kappa(1,t)];
        cd jpg
        eval(['saveas(H,''',filenam,'_',num2str(segfirst),'_',num2str(seglast),'_ZSHIFT.jpg'')']);
        cd ..
        fprintf(fid,'%s %f %f %f %f %f %f \n',...
            filenam,segfirst,seglast,...
            zmaxloc,zmax,hang,hmvl);      % write text file
    end     % end of "segment cycle"
    if ~isempty(ZMaxLoc(find(~isnan(ZMaxLoc))))
        ZMean(o) = b_mean_nonnan(ZMaxLoc);
        ZStd(o) = b_std_nonnan(ZMaxLoc);
    else
        ZMean(o) = NaN;
        ZStd(o) = NaN;
    end
    
% Save
    cd mat
    str = [filenam '_ZSHIFT'];
    save(str,'ZMaxLoc','ZMax','Ang','Mvl')
    str = [filenam '_ANG'];
    save(str,'Ang','Mvl')
    cd ..
    
    H2 = figure;
    plot(ZMaxLoc/sr,'.')
    cd jpg
    eval(['saveas(H2,''',filenam,'_WHOLEZ.fig'')']);
    cd ..
    close(H2)
    
    waitbar(o/sf)
end

AllMean = b_mean_nonnan(ZMean);
AllStd = b_std_nonnan(ZMean);
cd mat
str = 'all_ZSHIFT';
save(str,'ZMean','ZStd','AllMean','AllStd')

% Correlation between phase angles and 'zmaxloc'
fn1 = find(isnan(CumAng));
CumAng(fn1) = [];
CumZMaxLoc(fn1) = [];
CumKappa(fn1) = [];
fn2 = find(isnan(CumZMaxLoc));
CumAng(fn2) = [];
CumZMaxLoc(fn2) = [];
CumKappa(fn2) = [];
[b,bint,r,rint,stats] = regress(CumZMaxLoc',[ones(length(CumZMaxLoc),1) sin(CumAng') cos(CumAng')]);
AngZReg = stats;    % circular-linear regression
[b,bint,r,rint,stats] = regress(CumZMaxLoc',[ones(length(CumZMaxLoc),1),CumKappa']);
KappaZReg = stats;

save RegData AngZReg KappaZReg AngZReg CumAng
cd ..

% Finalize
datum = datevec(datestr(now));
save list list datum
cd(mmm)
fclose(fid);
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
function [zmaxloc,zmax,hang,hmvl] = zshift(H,eeg,vdisc,filenam,titlestr,segfirst,seglast,sr)

% Filtering EEG
nqf = sr / 2;
flt = fir1(4096,5/nqf,'low');      % lowpass filtering on 5 Hz
feeg = filtfilt(flt,1,eeg);

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
Z = [];
p = [];
T = [-1*sr:10:1*sr];
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
plot(T/sr,Z)

sl = 0.005;  % level of significance
siglev = -1 * log(sl);
hold on
line(xlim,[siglev siglev],'Color','k','LineStyle','--');

zmax = max(Z);
zmxlc = find(Z==zmax);
x_lim = xlim;
y_lim = ylim;
try
    text(x_lim(1)+2*(x_lim(2)-x_lim(1))/3,y_lim(1)+2*(y_lim(2)-y_lim(1))/3,num2str(T(zmxlc(1))));
end
tt = [titlestr ' ' num2str(segfirst) ' ' num2str(seglast)];
title(tt)
hold off

% Skip file, if phase connection is not significant at the maximum localization
if Z(zmxlc) < siglev
    zmaxloc = NaN;
else
    zmaxloc = T(zmxlc);
end