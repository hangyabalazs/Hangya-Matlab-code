function st = b_backup(last_backup)
%BACKUP   Creates backup files for Matlab work.
%   STATUS = BACKUP(LAST_BACKUP) saves all files in work directory that
%   have been modified since last backup to a destination directory
%   specified in the program code. ST = 1 if all files were successfully
%   copied and 0 otherwise.
%
%   See also FINISH.

% Creating machine independent work path
% work_path = fullfile(matlabroot,'work\Balazs\');
work_path = 'D:\Dropbox\code\MATLAB_work\';

% Destination
% global DATAPATH
% destroot = [DATAPATH 'Matlab work backup\Balazs\'];
% destroot = '\\science\Kepecs\Balazs\Matlab_backup\';
destroot = 'h:\Matlab_AUTOBACKUP2\';

% Copy
st = 1;   % status variable
dir_list{1} = work_path;
while ~isempty(dir_list)
    files = dir(dir_list{1});
    files = files(3:end);
    lf = length(files);
    for i = 1:lf
        if ~files(i).isdir
            dv = datevec(files(i).date);
            df = dv - last_backup;
            df(end+1) = -1;
            fdf = find(df);
            if df(fdf(1)) > 0
                source = fullfile(dir_list{1},files(i).name);
                destspec = dir_list{1}(length(work_path)+1:end);
                dest = fullfile(destroot,destspec);
                cmnd = ['copyfile(''',source,''',''',dest,''');'];
                try
                    eval(cmnd);
                catch
                    st = 0;
                    disp(['Unable to copy: ' source]);
                end
            end
        else
            %put new directory on the directory list
            new_dir = fullfile(dir_list{1},files(i).name);
            dirspec = dir_list{1}(length(work_path)+1:end);
            ff = fullfile(destroot,dirspec);
            cmnd = ['status = mkdir(''',ff,''',''',files(i).name,''');'];
            eval(cmnd);
            if isempty(find(strcmp(dir_list,new_dir)))
                dir_list{end+1} = new_dir;
            end
        end
    end
    dir_list(1) = [];
end