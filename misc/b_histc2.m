function nn = b_histc(y,bins)
%HISTC2   Histogram count.
%   HISTC2 does exactly the same as HISTC using another method. It is
%   usually slower than HISTC, but uses less memory.
%
%   See also HIST, HIST2 and HISTC.

lenb = length(bins);

% Input argument check
if isempty(bins) | lenb == 1
    error('Input argument ''bins'' has to contain at least 2 elements.')
end

dd = diff(bins);
df = find(dd==0);
if ~isempty(df)
    error('12')
end

% Histogram count
nn = zeros(1,lenb-1);
nn(1) = length(find(y<bins(2)));
for k = 2:lenb-2
    nn(k) = length(find(bins(k)<=y&y<bins(k+1)));
    k
end
nn(lenb-1) = length(find(bins(lenb-1)<=y));

leny = size(y,1) * size(y,2);
if sum(nn) ~= leny
    error('22')
end