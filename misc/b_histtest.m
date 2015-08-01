function b_histtest
%HISTTEST   Tests memory conserving version of histc.
%   HISTTEST tests HISTC2.
%
%   See also HISTC2 and HIST2.

for t = 1:1000
    k1 = fix(rand(1,1)*100) + 1;
    k2 = fix(rand(1,1)*15);
    k2 = max(k2,2);
    y = fix(rand(1,k1)*100);
    bins = zeros(1,k2);
    bins(1) = min(y);
    bins(end) = max(y);
    b = rand(1,length(bins)-2);
    b = b * (bins(end)-bins(1)) + bins(1);
    bins(2:end-1) = sort(fix(b));
    bins = dropequals(bins);
    if isempty(bins) | length(bins) == 1
        continue
    end
    nn1 = b_histc(y,bins);
    nn2 = histc(y',[-inf bins],1);
    nn2(2,:) = nn2(2,:) + nn2(1,:);
    nn2(end-1,:) = nn2(end-1,:) + nn2(end,:);
    nn2 = nn2(2:end-1,:)';
    if ~isequal(nn1,nn2)
        error('HISTTEST FAILED!');
    end
    if rem(t,100) == 0
        t
    end
end

% --------------------------------------------------------------------------------------------
function out = dropequals(bins)
binsold = bins;
for i = 1:length(bins)
    fnd = find(bins==binsold(i));
    lfnd = length(fnd);
    for j = lfnd:-1:2
        bins(fnd(j)) = [];
    end
end
out = bins;