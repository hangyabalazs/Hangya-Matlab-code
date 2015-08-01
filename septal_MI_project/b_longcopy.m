function b_longcopy(long,source)
%LONGCOPY   Copies long raw data files to another folder.
%   LONGCOPY(LONG) copies the files specified in LONG to the destination directory
%   which you can modify in the program code. LONG is usually created by THETASELECTORRUN
%   and contains the name of the files that were too long for THETASELECTORRUN.
%   You can analyse these files using THETASELECTORRUN_LONG.
%
%   LONGCOPY(LONG,SOURCE) uses SOURCE as source directory.
%
%   See also THETASELECTORRUN and THETASELECTORRUN_LONG.

% Input arguments check
error(nargchk(1,2,nargin));

% Searching for the matching file
global DATADIR
global DATAPATH
if nargin < 2
    source = DATADIR;
end
dest = [DATAPATH,'DATA\analysenow3\'];
k = dir(source);
j = zeros(length(k),1);
cell_j = cell(1,length(k));
for m = 3:length(k),
    cell_j{m} = k(m).name(1:6);
end

% Copy
for i = 1:length(long)
    j = strcmp(cell_j,long{i});
    if isempty(find(j)) == 0,
        filename = k(find(j)).name;
        pont = findstr(filename,'.');
        filenam = filename(1:pont(1)-1);
    end
    fn = fullfile(source,filename);
    eval(['copyfile(''',fn,''',''',dest,''');']);
end