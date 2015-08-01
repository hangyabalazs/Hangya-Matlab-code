function [files2, files2_short] = b_filelist(inpdir)
%FILELIST    List of files in directory.
%   FILES = FILELIST(INPDIR) returns the list of files in INPDIR directory.
%   [FILES, SHORT] = FILELIST(INPDIR) returns the first six characters of
%   filenames in SHORT.
%
%   See also DIR.

% Input argiment check
error(nargchk(1,1,nargin))

% List of filenames
files = dir(inpdir);
files = files(3:end);
files2 = struct([]);
files2_short = {};
for i = 1:length(files)
    if ~files(i).isdir
        if isempty(files2)
            files2 = files(i);
        else
            files2(end+1) = files(i);
        end
        files2_short{end+1} = files(i).name(1:min(6,length(files(i).name)));
    end
end