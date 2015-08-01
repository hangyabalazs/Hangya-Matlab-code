function noseg
%NOSEG   Analysis of theta and non-theta segments (NO project).
%   NOSEG calculates the frequency of non-theta segments in the last 1800 s
%   for pre-injection files or in the first 3600 s in post-injection files.
%   The results are saved in an Excel file.
%
%   See also NOSEG2.

% Directories
global DATADIR
global DATAPATH
inpdir = [DATADIR 'NO_matfiles\'];
inpdir2 = [DATAPATH 'NO\Wavelet\Segments\'];
resdir = [DATAPATH 'NO\Seg\'];
tblfile = [DATAPATH 'NO\seg_data.xls'];
% mm = pwd;
% cd(resdir)

% Filelist
[files files_short] = filelist(inpdir);
sf = length(files_short);

% Progress indicator
wb = waitbar(0,'Running NOSEG...','Position',[360 250 275 50]);
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
ntfrs = zeros(1,sf);
nls = zeros(1,sf);
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
    
    ntf = ThetaSegments(1,:);   % frequency of non-theta segments
    ntl = ThetaSegments(2,:);
    ntlen = (ntl - ntf) / sr;
    if strcmp(id2,'pre')
        nt = ntf(ntf>len-1800*sr);
        ntfr = length(nt) / 1800;
        nl = ntlen(ntf>len-1800*sr);
    elseif strcmp(id2,'post')
        nt = ntf(ntl<1800*sr);
        ntfr = length(nt) / 1800;
        nl = ntlen(ntl<1800*sr);
    end
    ntfrs(o) = ntfr;
    nls(o) = sum(nl) / 1800;
end

% Save
xlswrite([resdir 'nonthetafreq5.xls'],names','sheet1','A1')
xlswrite([resdir 'nonthetafreq5.xls'],nls','sheet1','B1')

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