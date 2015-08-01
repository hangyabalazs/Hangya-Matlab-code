function b_rescaleaxis(axis,scale)
%RESCALEAXIS   Changes ticklabeles to 'frequency'.
%   RESCALEAXIS(AXIS,SCALE) requires input parameter SCALE containing the correspondig
%   frequency value for each index. It changes the ticklabels of AXIS to the
%   frequency values.
%
%   See also RETIMEFIG and RESCALEFIG.

%Input arguments check
error(nargchk(2,2,nargin));

% Get label
if strcmp(axis,'X')
    a0 = get(gca,'XTick')';
elseif strcmp(axis,'Y')
    a0 = get(gca,'YTick')';
else
    error('First input argument must be ''X'' or ''Y''.')
end

% Interpolation
lo0 = floor(a0);
if ~isequal(a0,lo0);
    hi0 = ceil(a0);
    iv = hi0 - lo0;
    liv = a0 - lo0;
    fiv = find(iv);
    r = zeros(size(iv));
    r(fiv) = liv(fiv) ./ iv(fiv);
    milo = lo0 < 1 | lo0 > length(scale);
    mihi = hi0 < 1 | hi0 > length(scale);
    lo1 = zeros(size(lo0));
    hi1 = zeros(size(hi0));
    lo1(find(milo)) = NaN;
    hi1(find(mihi)) = NaN;
    lo1(find(~milo)) = scale(lo0(find(~milo)));
    hi1(find(~mihi)) = scale(hi0(find(~mihi)));
    a1 = lo1 + r .* (hi1 - lo1);
else
    mia = a0 < 1 | a0 > length(scale);
    a1 = zeros(size(a0));
    a1(find(mia)) = NaN;
    a1(find(~mia)) = scale(a0(find(~mia)));
end
a1 = round(a1*100) / 100;

% Set label
na1 = num2str(a1);
[m n] = size(na1);
rna = reshape(na1',1,m*n);
ni = findstr(rna,'NaN');
rna(ni) = ' ';
rna(ni+1) = ' ';
rna(ni+2) = ' ';
na1 = reshape(rna,n,m)';
if strcmp(axis,'X')
    set(gca,'XTickLabel',na1)
elseif strcmp(axis,'Y')
    set(gca,'YTickLabel',na1)
end