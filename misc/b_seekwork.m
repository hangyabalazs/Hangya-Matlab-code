function fullname = b_seekwork(fcn_name)
%SEEKWORK   Creates full work dir. filename from partial file name.
%   FULLNAME = SEEKWORK(FCN_NAME) looks for valid function names in the
%   work directory and its subdirectories.
%
%   See also WHICH.

%Input argument check.
error(nargchk(1,1,nargin));
if ~isstr(fcn_name)
    error('Input must be strings.');
end

%Create machine independent work path
work_path = fullfile(matlabroot,'work');
lw = length(work_path);

%Searching...
[pathname, fname, extension] = fileparts(fcn_name);
if isempty(extension)
    if exist(pathname,'dir')   
        %have full path and name but no extension
        if strncmpi(pathname,work_path,lw) == 0     %check if path is work dir
            error('File must be in the work directory.');
        end
        if exist([fcn_name '.m'])       %look for .m, .p and mex files
            extension = '.m';
        elseif exist([fcn_name mexext])
            extension = mexext;
        elseif exist([fcn_name '.p'])
            extension = 'p';
        else
            error(['File ' fcn_name ' cannot be located or is not an ' ...
                    'M/MEX/P file.']);
        end
        fullname = [fullfile(pathname,fname) extension];
    else
        %have no path and no extension
        fullname = which(fcn_name,'-all');      %search
        if isempty(fullname)        %no result
            error(['Cannot locate file ' fcn_name '.']);
        end
        fnd = find(strncmpi(fullname,work_path,lw));
        if isempty(fnd)     %check if any path is work dir
            error('Input file must be in the work directory.');
        elseif length(fnd) > 1      %decision if more results
            imp = cell(length(fnd),1);
            for r = 1:length(fnd)
                imp{r} = [num2str(r) '   ' fullname{fnd(r)}];
            end
            char_imp = char(imp);
            disp(['Which one do you want to examin?']);
            disp(char_imp);
            npt = input(['Choose a number from 1 to ' num2str(length(fnd))...
                    '! ']);
            if isempty(find([1:length(fnd)]==npt))
                error('Unexpected answer.');
            end
            fullname = fullname{fnd(npt)};
        end
    end
elseif exist(pathname,'dir')
    %got an extension and a full path
    if strncmpi(pathname,work_path,lw) == 0     %check if path is work dir
        error('File must be in the work directory.');
    end
    fullname = fcn_name;
elseif any(strcmp(extension,{'.m' '.p' ['.' mexext]}))      %check the extension
    %got an extension but no path
    fullname = which(fcn_name,'-all');      %search
    if isempty(fullname)        %no result
        error(['Cannot locate file ' fcn_name '.']);
    end
    fnd = find(strncmpi(fullname,work_path,lw));
    if isempty(fnd)     %check if path is work dir
        error('Input file must be in the work directory.');
    elseif length(fnd) > 1      %decision if more results
        imp = cell(length(fnd),1);
        for r = 1:length(fnd)
            imp{r} = [num2str(r) '   ' fullname{fnd(r)}];
        end
        char_imp = char(imp);
        disp(['Which one do you want to examin?']);
        disp(char_imp);
        npt = input(['Choose a number from 1 to ' num2str(length(fnd))...
                '! ']);
        if isempty(find([1:length(fnd)]==npt))
            error('Unexpected answer.');
        end
        fullname = fullname{fnd(npt)};
    end
else
    error(['File ' fcn_name ' is not an allowed file type (M/P).']);  
end
%check if result is valid
if iscell(fullname)
    fullname = fullname{1};
end
[pathname fcn_name extension] = fileparts(fullname);
if ~any(strcmp(extension,{'.m' '.p' ['.' mexext]})) | exist(fullname) == 0
    error(['File ' fcn_name ' cannot be located or is not an ' ...
            'M/MEX/P file.']);
end

%   List of files calling SEEKWORK:
%   b_fundep