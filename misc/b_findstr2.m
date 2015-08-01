function fs = b_findstr2(str1,varargin)
%FINDSTR2   Find a cell array of strings within a single string.
%   FS = FINDSTR2(STR1, STR2, ...) returns indeces of startingpoints of
%   STR2 or other optional input strings within STR1. For example,
%   findstr2('Life is hard if you dont use this program.','if','ha')
%   returns [2, 9, 14].
%
%   See also FINDSTR and STRREAD2.

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