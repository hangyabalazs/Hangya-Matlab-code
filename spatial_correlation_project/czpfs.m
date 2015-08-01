function R = czpfs(rm1,rm2)
%CZPFS   Place Field Similarity.
%   R = CZPFS(RM1,RM2) calculates linear correlation between rate maps RM1
%   and RM2.
%
%   See also CZPLACEANALYSIS and CZPFS_MOD.

% Find non-NaNs
r1 = rm1(~isnan(rm1)&~isnan(rm2));
r2 = rm2(~isnan(rm1)&~isnan(rm2));

% Linear correlation
corrmtx = corrcoef(r1,r2);
R = corrmtx(2);