function [list,error_files] = b_fastfundep(fcn_name)
%FASTFUNDEP Finds work dir. functions that call a specific work dir. function.
%   [LIST,ERROR_FILES] = FASTFUNDEP(FCN_NAME) creates a LIST of functions that 
%   are in the work directory or in its subdirectory and call FCN_NAME. ERROR_FILES
%   contains the list of files that produced an error during examination. These
%   problems can arise from failure of opening the file with FOPEN.
%
%   FASTFUNDEP, despite of its name, is slower then FUNDEP and is less detailed,
%   but doesn't need FUNCTIONSCALLED and temporary file.
%
%   Note, that FCN_NAME should be in the work directory or in its subdirectory!
%
%   You should give the input function name without path or extension!
%
%   See also DEPFUN, PROFILE and FUNDEP.

%Input argument check
error(nargchk(1,1,nargin));
if ~isstr(fcn_name)
    error('Input must be strings.');
end

%Creating machine independent work path
work_path = fullfile(matlabroot,'work\Balazs');
disp(['Input directory: ' work_path])
lw = length(work_path);

%Create full work dir. filename of input function
fullname = [work_path '\' fcn_name '.m'];

%Initializing some variables
dir_list{1} = work_path;
list = {};
error_files = {};

%Seek for the calling functions in work directories
while ~isempty(dir_list)
    trace_list{1} = fullname;
    trace_list_junk = {};
    while ~isempty(trace_list)
        [pathname fname extension] = fileparts(trace_list{1});
        files = dir(dir_list{1});
        size_files = size(files,1);
        for i = 3:size_files
            if files(i).isdir == 0
                examin_fcn = fullfile(work_path,files(i).name);
                [examin_path examin_fname examin_ext] = fileparts(examin_fcn);
                %continue, if examined function is already on the list or is
                %equal to the input function (to avoid problem of self calling)
                if ~isempty(find(strcmp(list,examin_fcn))) | strcmp(fname,examin_fname)
                    continue
                end
                if any(strcmp(examin_ext,{'.m' '.p' ['.' mexext]}))
                    %check out code looking for function calls
                    try
                        result = b_fastfunctionscalled(fname,examin_fcn);
                        %if there is a result, put it on the list and the trace list
                        if result
                            if isempty(find(strcmp(trace_list,examin_fcn)))...
                                    & isempty(find(strcmp(trace_list_junk,examin_fcn)))
                                trace_list{end+1} = examin_fcn;
                                list{end+1} = examin_fcn;
                            end
                        end
                    catch
                        %if there was an error while running FASTFUNCTIONSCALLED, 
                        %put the file name on error files list
                        if isempty(find(strcmp(error_files,examin_fcn)))
                            error_files{end+1} = examin_fcn;
                        end
                    end
                end
            else
                %put new directory on the directory list
                new_dir = fullfile(work_path,files(i).name);
                if isempty(find(strcmp(dir_list,new_dir)))
                    dir_list{end+1} = new_dir;
                end
            end    
        end
        trace_list_junk{end+1} = trace_list{1};
        trace_list(1) = [];
    end
    dir_list(1) = [];
end

%Display the list of error files
if ~isempty(error_files)
    disp(' ');
    disp('These functions could not be parsed: ');
    for r = 1:length(error_files)
        disp(error_files{r});
    end
end

%Transpose list and error_files
list = list';
error_files = error_files';