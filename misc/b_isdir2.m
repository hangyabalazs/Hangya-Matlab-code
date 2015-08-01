function logic = b_isdir2(dirname)
%ISDIR2   True for subdirectories of current directory.
%   ISDIR(DIR) returns a 1 if DIR is a subdirectory in the current directory
%   and 0 otherwise.
%
%   See also EXIST and ISDIR.

% Input argument check
error(nargchk(1,1,nargin));

% Output
d = dir;
lend = length(d);
n = cell(1,lend);
for i = 1:lend
    n{i} = d(i).name;
end
logic = ~isempty(find(strcmp(dirname,n)));