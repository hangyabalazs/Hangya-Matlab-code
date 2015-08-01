function data = b_load_data(ff)
%LOAD_DATA   Load raw data.
%   DATA = LOAD_DATA(FF) loads row data with the name stored in FF (string input)
%   and assigns DATA the loaded matrix. It also makes some necessary transformations.
%
%   See also LOAD, SUBSREF and IN3.

data = load(ff);   % load data
if isstruct(data)
    field = fieldnames(data);
    s = struct('type','.','subs',field);
    data = subsref(data,s);
end
if size(data,2) == 1
    data = [data(1:length(data)/2,1) data((length(data)/2)+1:length(data),1)];
end