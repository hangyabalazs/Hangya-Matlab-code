function b_nontheta(f,newstep,wdir,where2)
%NONTHETA   Non-theta segment selector.
%   NONTHETA(F,NEWSTEP) loads THETASELECTORRUN output matrices (OM) and saves 
%   'NonThetaSegments' matrix and 'SharpWaveSegments' matrix, who are 2-N matrices
%   containing first points of non-theta or sharpwave segments in their first row and
%   last points of non-theta or sharpwave segments in their second row. F should be the
%   scalevector used in WAVELET and NEWSTEP should be the constant which characterizes 
%   downsampling of the EEG (sampled on 10 kHz).
%
%   NONTHETA uses four directories: two data directories (thetaselection data directory
%   and raw data directory) and two result directories (nontheta segments' directory and
%   sharpwave segments' directory). You are able to modify them through editing the program
%   code, or you can use NONTHETA(F,NEWSTEP,WDIR) syntax, where WDIR is the working directory
%   of NONTHETA function (in which the result directories are defined automatically). WDIR
%   should contain thetaselection data directory - usually THETASELECTORRUN does this job.
%   You can add raw data directory as input argument with NONTHETA(F,NEWSTEP,WDIR,RAWDIR)
%   syntax.
%
%   Theta band: 2.5 - 6 Hz.
%
%   Note: not longer than 6144 segments are skipped (unable to filter)!
%
%   See also THETASELECTORRUN, NONTHETA and THETA_LONG.

% Input arguments check
error(nargchk(2,4,nargin));

% Define directories
global DATAPATH
global DATADIR
if nargin < 3
    wdir = [DATAPATH,'HCN\Wavelet2\'];
end
where = [wdir 'thetaselection\matrix\'];    %Thetaselection data directory
create_subdir(wdir)   % create subdirectories
noth_dir = [wdir 'nontheta_segments\'];  %Nontheta segments' directory
sharp_dir = [wdir 'sharpwave_segments\'];  %Sharpwave segments' directory
if nargin < 4
    where2 = ['f:\raw_data\hcn\all\'];    %Raw data directory
end

% Import
files = dir(where);
files = files(3:end);
sf = length(files);
mm = pwd;

files2 = dir(where2);
ns_long = cell(1,length(files2));
ns_short = cell(1,length(files2));
for i = 3:length(files2)
    ns_long{i} = files2(i).name;
    fs = findstr(ns_long{i},'_');
    ns_short{i} = ns_long{i}(1:fs(2)-1);
end

wb = waitbar(0,'Running NONTHETA...','Position',[360 250 275 50]);    %Progress indicator
global WB
WB(end+1) = wb;

% Loading THETASELECTORRUN output matrix
for o = 1:sf
    filenam = files(o).name;
    ff = [where filenam];
    load(ff)
    maxloc = Out(1,:);
    
% Loading data
    fs = findstr(filenam,'_');
    fsp = findstr(filenam,'.');
    scp = strcmp(filenam(fs(2)+1:fsp(1)-1),ns_short);
    fnd = find(scp);
    fln = ns_long{fnd};
    ff2 = [where2 fln];
    eeg = [];
    load(ff2);
    if isempty(eeg)
        eeg = data(:,1)';
    end
    
% Finding complement of theta segments
    fnd = find(f>6);
    pwind1 = fnd(end);
    fnd = find(f<2.5);
    pwind2 = fnd(1);
    fml = find((maxloc<=pwind1)|(maxloc>=pwind2));
    
% Finding non-theta segments
    if ~isempty(fml)
        thcomp = segm(fml,newstep);     %Theta complementer set
        [noth,shwa] = segm2(thcomp,eeg);       %Non-theta segments, sharp waves
    end
    
% Save
    NonThetaSegments = noth;
    SharpWaveSegments = shwa;
    noth = [];
    shwa = [];
    fs = findstr(filenam,'_');
    fln1 = filenam((fs(2)+1):(fs(3)-1));
    fln2 = filenam((fs(3)+1):end-4);
    cd(noth_dir)
    eval(['save(''NONTHETA_SEGMENTS_',fln1,'_',fln2,'.mat'',''NonThetaSegments'')']);
    cd ..
    cd(sharp_dir)
    eval(['save(''SHARPWAVE_SEGMENTS_',fln1,'_',fln2,'.mat'',''SharpWaveSegments'')']);
    cd ..
    
    waitbar(o/sf)   %Progress indicator
end
close(wb);   %Close progress indicator

% ----------------------------------------------------------------------------------------
function thcomp = segm(ip,newstep)
drml = diff(ip);
fdr = find(drml>1);
lenfdr = length(fdr);
preth = zeros(2,lenfdr+1);
preth(1,1) = ip(1);
for t = 1:lenfdr
    preth(2,t) = ip(fdr(t));
    preth(1,t+1) = ip(fdr(t)+1);
end
preth(2,end) = ip(end);
thcomp = preth .* newstep;

thc1 = thcomp(1,:);
thc2 = thcomp(2,:);
difthc = thc2 - thc1;
fd = find(difthc<=6144);         %leaving segments not longer than 0,6144 sec.
thc1(fd) = [];
thc2(fd) = [];
thcomp = [thc1; thc2];

% ----------------------------------------------------------------------------------------
function [noth, shwa] = segm2(thcomp,eeg)
noth1 = [];
noth2 = [];
shwa1 = [];
shwa2 = [];
for i = 1:size(thcomp,2)
    noth1(end+1) = thcomp(1,i);
    ifrst = max(1,thcomp(1,i));
    iscnd = min(thcomp(2,i),length(eeg));
    [bg fin] = b_sspw_nontheta(eeg(ifrst:iscnd));
    shwa1 = [shwa1 thcomp(1,i) + bg];
    shwa2 = [shwa2 thcomp(1,i) + fin];
    if ~isempty(bg)
        for j = 1:length(bg)+1
            if j == 1
                noth2(end+1) = thcomp(1,i) + bg(1);
            elseif j == length(bg) + 1
                noth1(end+1) = thcomp(1,i) + fin(end);
            else
                noth1(end+1) = thcomp(1,i) + fin(j-1);
                noth2(end+1) = thcomp(1,i) + bg(j);
            end
        end
    end
    noth2(end+1) = thcomp(2,i);
end
noth = [noth1; noth2];
shwa = [shwa1; shwa2];

% ----------------------------------------------------------------------------------------
function create_subdir(wdir)

mm = pwd;
cd(wdir)
if ~b_isdir2('nontheta_segments')
    mkdir nontheta_segments
end
if ~b_isdir2('sharpwave_segments')
    mkdir sharpwave_segments
end
cd(mm)