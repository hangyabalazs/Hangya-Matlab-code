function b_zshiftrun
%ZSHIFTRUN    Runs ZSHIFT on a set of files.
%   ZSHIFTRUN uses an input and a result directory: user can set them
%   through editing the program code.
%
%   Maximum localizations of Z arrays for theta and non-theta segments,
%   there means and standard deviations for one cell and along all cells
%   are getting saved in mat files. The plots are getting saved in jpg
%   files.
%
%   ZSHIFTRUN uses its own ZSHIFT instead of ZSHIFT2.
%
%   See also ZSHIFT2.

% Input argument check
error(nargchk(0,0,nargin))

% Directories
global DATAPATH
global DATADIR
where = [DATADIR,'all/'];    %Here are the discriminated data files
% where = [DATAPATH,'DATA/analysenow5/'];    %Here are the discriminated data files
files = dir(where);
files = files(3:end);
sf = length(files);
mmm = pwd;
cd([DATAPATH,'Zshift/']);  %Here are the results
create_subdir;         %create subdirectories

% Progress indicator
wb = waitbar(0,'Running ZSHIFTRUN...','Position',[360 250 315 50]);    %Progress indicator
global WB
WB(end+1) = wb;

% Open figure
H = figure;

% Import
list = {};      %initialize list of analysed files
ThetaMean = zeros(1,sf);
ThetaStd = zeros(1,sf);
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
    ff = fullfile(DATAPATH,'Wavelet2/theta_segments',fln);
    load(ff)
    
    segments = ThetaSegments;
    if ~isempty(segments)
        segments = short_killer(segments);      % leave 200 - 200 ms edges then skip segments
    end

% Run 'zshift' on theta intervals
    th_index = size(segments,2);
    ThetaZMaxLoc = zeros(1,th_index);
    for t = 1:th_index      % "THETA SEGMENT CYCLE"
        segfirst = segments(1,t);
        seglast = segments(2,t);
        eeg_input = eeg(segfirst:seglast);
        vdisc_input = vdisc(find(vdisc>=segfirst&vdisc<=seglast)) - segfirst + 1;
        
        zmaxloc = zshift(H,eeg_input,vdisc_input);
        ThetaZMaxLoc(t) = zmaxloc;
        cd jpg
        eval(['saveas(H,''THETA_',filenam,'_',num2str(segfirst),'_',num2str(seglast),'_ZSHIFT.jpg'')']);
        cd ..
        if zmaxloc == NaN   % add skipedd file to a list
            
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
    ff = fullfile(DATAPATH,'Wavelet2/nontheta_segments',fln);
    load(ff)
    
    segments = NonThetaSegments;
    if ~isempty(segments)
        segments = short_killer(segments);      % leave 200 - 200 ms edges then skip segments
    end

% Run 'zshift' on non-theta intervals
    no_index = size(segments,2);
    NonThetaZMaxLoc = zeros(1,no_index);
    for t = 1:no_index      % "THETA SEGMENT CYCLE"
        segfirst = segments(1,t);
        seglast = segments(2,t);
        eeg_input = eeg(segfirst:seglast);
        vdisc_input = vdisc(find(vdisc>=segfirst&vdisc<=seglast)) - segfirst + 1;
        
        zmaxloc = zshift(H,eeg_input,vdisc_input);
        NonThetaZMaxLoc(t) = zmaxloc;
        cd jpg
        eval(['saveas(H,''NONTHETA_',filenam,'_',num2str(segfirst),'_',num2str(seglast),'_ZSHIFT.jpg'')']);
        cd ..
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
    save(str,'ThetaZMaxLoc','NonThetaZMaxLoc')
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
cd ..

% Finalize
datum = datevec(datestr(now));
save list list datum
cd(mmm)
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
function zmaxloc = zshift(H,eeg,vdisc)

% Hilbert transformation of the eeg
ahee = angle(hilbert(eeg));
aheedeg = ahee * (180 / pi);

% Shift unit
Z = [];
p = [];
T = [-10000:10:10000];
for t = T
    vt = vdisc(find(vdisc+t>0&vdisc+t<length(eeg)))+t;  % shift unit
    if isempty(vt)      % Skip segment, if all spikes fall out of range after unit shift
        zmaxloc = NaN;
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

zmxlc = find(Z==max(Z));
x_lim = xlim;
y_lim = ylim;
text(x_lim(1)+2*(x_lim(2)-x_lim(1))/3,y_lim(1)+2*(y_lim(2)-y_lim(1))/3,num2str(zmxlc(1)));
hold off

% Skip file, if phase connection is not significant at the maximum localization
if Z(zmxlc) < siglev
    zmaxloc = NaN;
else
    zmaxloc = T(zmxlc);
end