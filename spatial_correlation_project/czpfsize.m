function pfsz = czpfsize(rm)
%CZPFSIZE   Place field size.
%   PS = CZPFSIZE(RM) calculates place field size for a rate map (RM) by
%   calculating the maximal area (in pixels) where firing rate is greater 
%   than mean + sd.
%
%   See also CZPLACEANALYSIS.

% Threshold criterium (mean + sd)
% ri = rm .* (zero2nan(double(rm>b_mean_nonnan(rm)+2*b_std_nonnan(rm))));
ri = rm .* (zero2nan(double(rm>b_max_nonnan(rm)*0.7)));

% Find connected components
L = bwlabel(nan2zero(ri));
nc = max(L(:));
sc = zeros(1,nc);
for k = 1:nc;
    sc(k) = length(find(L==k));
end
pfsz = max(sc);