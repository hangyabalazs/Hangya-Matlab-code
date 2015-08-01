function [list,eval_strings,error_files] = b_fundep(fcn_name)
%FUNDEP Finds work dir. functions that call a specific work dir. function.
%   [LIST,EVAL_STRING,ERROR_FILES] = FUNDEP(FCN_NAME) creates a LIST of 
%   functions that are in the work directory or in its subdirectory and
%   call FCN_NAME. EVAL_STRING returns lines that contain EVAL, FEVAL,
%   EVALIN functions. ERROR_FILES contains the list of files that produced
%   an error while being parsed. Parsing problems can arise from MATLAB 
%   syntax errors.
%
%   Note, that FCN_NAME should be in the work directory or in its subdirectory!
%
%   See also DEPFUN, PROFILE and FASTFUNDEP.

% Input argument check
error(nargchk(1,1,nargin));
if ~isstr(fcn_name)
    error('Input must be strings.');
end

% Creating machine independent work path
work_path = fullfile(matlabroot,'work');
lw = length(work_path);

% Make sure tempdir is on the path
if isempty(findstr(tempdir,matlabpath))
  addpath(tempdir)
  added_temp = 1;
else
  added_temp = 0;
end
% Generate a temp m-file name that is used to hold evaluated strings
% for compilation.  Test temp file to verify it can be opened in write
% mode and deleted.  Unfortunately the delete command does not provide
% an error return value so I'm catching the warning and then erroring out.
lastwarn('');
tmpname = [tempname '.m'];
fid = fopen(tmpname, 'w');
fclose(fid);
if exist(tmpname), delete(tmpname); end
if strcmp(lastwarn,'File not found or permission denied')
    error(['Unable to run DEPFUN due to inability to write to and remove ' ...
            'temporary file (' tmpname ').']);
end

% Create full work dir. filename of input function
fullname = b_seekwork(fcn_name)

% Initializing some variables
dir_list{1} = work_path;
list = {};
eval_strings = [];
error_files = {};

cbparam = {};
HGVariables = {};
HGVariableStrings = {};

% Seek for the calling functions in work directories
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
                if ~isempty(find(strcmp(list,examin_fcn)))| strcmp(fname,examin_fname)
                    continue
                end
                if any(strcmp(examin_ext,{'.m' '.p' ['.' mexext]}))
                    %check out code looking for function calls - using the builtin
                    %function FUNCTIONSCALLED
                    try
                        [names, dyn_syms, problem_lines, java_classes, java_imports,...
                                java_maybe_classes] = functionscalled(examin_fcn, tmpname,...
                            cbparam, HGVariables, HGVariableStrings);
                        if exist(tmpname)
                            delete(tmpname);
                        end
                        exist_calls = strcmpi(names,fname);
                        fnd_calls = find(exist_calls);
                        %if there is a result, put it on the list and the trace list
                        if ~isempty(fnd_calls)
                            if isempty(find(strcmp(trace_list,examin_fcn)))...
                                    & isempty(find(strcmp(trace_list_junk,examin_fcn)))
                                trace_list{end+1} = examin_fcn;
                                list{end+1} = examin_fcn;
                            end
                        end
                        %lines containing evaluating functions
                        if ~isempty(problem_lines)
                            eval_strings(end+1).fcn_name = examin_fcn;
                            eval_strings(end).lineno = problem_lines;
                        end
                    catch
                        %if there was an error while running FUNCTIONSCALLED, put the
                        %file name on error files list
                        if exist(tmpname)
                            delete(tmpname);
                        end
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

% Be sure to clean up temp m-file
if exist(tmpname)
    delete(tmpname); 
end
if added_temp
    rmpath(tempdir)
end

% Display the list of error files
if ~isempty(error_files)
    disp(' ');
    disp('These functions could not be parsed: ');
    for r = 1:length(error_files)
        disp(error_files{r});
    end
end

% Transpose list and error_files
list = list';
error_files = error_files';