function noseg2
%NOSEG2   Analysis of theta and non-theta segments (NO project).
%   NOSEG2 calculates the amplitude (SD) of non-theta segments in the last
%   1800 s for pre-injection files or in the first 3600 s in post-injection
%   files. The results are saved in an Excel file.
%
%   See also NOSEG.

% Directories
global DATADIR
global DATAPATH
inpdir = [DATADIR 'NO_matfiles\'];
inpdir2 = [DATAPATH 'NO\Wavelet\Segments\'];
resdir = [DATAPATH 'NO\Seg\'];
tblfile = [DATAPATH 'NO\seg_data.xls'];

% Filelist
[files files_short] = filelist(inpdir);
sf = length(files_short);

% Progress indicator
wb = waitbar(0,'Running NOSEG2...','Position',[360 250 275 50]);
global WB
WB(end+1) = wb;

% Load segment data
[tbl0 tbl] = xlsread(tblfile);

% Main
sr = 5000;     % sampling rate
nqf = sr / 2;   % Nyquist frequency
newstep = 50;
dsr = sr / newstep;    % downsample on 100 Hz
names = cell(1,sf);
ntsds = zeros(1,sf);
for o = 1:sf
    fname = files(o).name;      % load data
    names{o} = fname;
    cmps = strread(fname(1:end-4),'%s','delimiter','_');
    titlestr = [cmps{1} ' ' cmps{2}];
    
    load([inpdir fname])
    eeg = hEEG.values;
    len = length(eeg);
    
    fn2 = [fname(1:end-4) '_SEGMENTS'];     % load segment identification
    load([inpdir2 fn2])
    inx = find(strcmp({tbl{:,1}},fname));
    id1 = tbl{inx,2};
    id2 = tbl{inx,3};
    
    ntf = NonThetaSegments(1,:);        % SD of non-theta segments
    ntl = NonThetaSegments(2,:);
    ntlen = (ntl - ntf) / sr;
    if strcmp(id2,'pre')
        ntinx = find(ntf>len-1800*sr);
    elseif strcmp(id2,'post')
        ntinx = find(ntl<1800*sr);
    end
    ntfi = ntf(ntinx);
    ntli = ntl(ntinx);
    lenin = length(ntinx);
    ntsd = zeros(1,lenin);
    for k = 1:lenin
        ntsd(k) = std(eeg(ntfi(k):ntli(k)));
    end
    ntsds(o) = mean(ntsd);
end

% Save
xlswrite([resdir 'nonthetasd.xls'],names','sheet1','A1')
xlswrite([resdir 'nonthetasd.xls'],ntsds','sheet1','B1')

% cd(mm)

% -------------------------------------------------------------------------
function [files2 files2_short] = filelist(inpdir)

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