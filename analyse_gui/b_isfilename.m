function L = b_isfilename(str)
%ISFILENAME   True if argument is a file name.
%   ISFILENAME(STR) returns a 1 if STR is a valid full file name (with path
%   and extension) and 0 otherwise.
%
%   See also EXIST, ISDIR and ISEMPTY.

[pat,nam,ext] = fileparts(str);
ne = [nam ext];
list = dir(pat);
lenl = length(list);
cell_list = cell(1,lenl-2);
for m = 3:lenl
    cell_list{m-2} = list(m).name;
end
ss = strcmp(ne,cell_list);
L = ~isempty(find(ss));