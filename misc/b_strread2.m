function cmps = b_strread2(str1,varargin)
%STRREAD2    Read  string using multiple delimiter.
%   CMPS = STREAD(STR1,STR2,...) reads the components of STR1 to CMPS cell
%   array using STR2,... as delimiters.
%
%   See also STRREAD and FINDSTR2.

% Input argument check
if nargin < 2
    error('At least 2 input arguments required.')
end

% Find string
R = [];
for i = 1:length(varargin);
    R = [R findstr(str1,varargin{i})];
end
fs = sort(R);
fs = [0 fs length(str1)+1];

% Read string
cmps = {};
for j = 1:length(fs)-1
    cmps{end+1} = str1(fs(j)+1:fs(j+1)-1);
end