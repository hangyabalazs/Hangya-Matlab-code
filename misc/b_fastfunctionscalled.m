function result = b_fastfunctionscalled(fun_called,fun_calls)
%FASTFUNCTIONSCALLED    Check if a work dir. function calls another.
%   RESULT = FASTFUNCTIONSCALLED(FUN_CALLED,FUN_CALLS) returns 1, if FUN_CALLES
%   calls FUN_CALLED, and 0 otherwise.
%
%   FASTFUNCTIONSCALLED, despite of its name, is slower then DEPFUN and is less detailed,
%   but doesn't need a temporary file.
%
%   Note, that the input functions should be in the work directory or in its
%   subdirectory!
%
%   You should give FUN_CALLED without path or extension!
%
%   See also FASTFUNDEP.

%Input argument check
error(nargchk(2,2,nargin));
if ~isstr(fun_called) | ~isstr(fun_calls)
    error('Input must be strings.');
end

%Create machine independent work path
work_path = fullfile(matlabroot,'work');

%Create full filename of fun_calls
[calls_path calls_fname calls_ext] = fileparts(fun_calls);
if isempty(calls_ext)
    if isempty(calls_path)
        fun_calls = [work_path '\' fun_calls '.m'];
    else
        fun_calls = [fun_calls '.m'];
    end
else
    if isempty(calls_path)
        fun_calls = [work_path '\' fun_calls];
    end
end

%Read the function
cell_lines = {};
fid = fopen(fun_calls,'r');
while ~feof(fid)
    line = fgetl(fid);
    if isempty(line) | strncmp(line,'%',1)...       %skip commented end empty
            | ~isempty(findstr('disp',line))...     %lines and lines containing
            | ~isempty(findstr('input',line))...    %some special statements
            | ~isempty(findstr('error',line))...
            | ~isempty(findstr('warning',line))...
            | ~isempty(findstr('eval',line))...
            | ~isempty(findstr('evalin',line))...
            | ~isempty(findstr('feval',line))
        continue
    end
    cell_lines{end+1} = line;
end

%Creating result using FINDSTR
if length(cell_lines) ~= 0
    for z = 1:length(cell_lines)
        result = ~(isempty(findstr(cell_lines{z},fun_called)));
        if result
            return
        end
    end
else
    result = 0;
end

%Another idea:
%FFF = textread('animation.m','%s','delimiter',...
%['\n',';','(',')',' ',',','/','*','+','-'],'whitespace','');

%Another idea:
%fid = fopen('c:\matlabR12\work\animation.m');
%F = fscanf(fid,'%c',inf);
%fclose(fid)