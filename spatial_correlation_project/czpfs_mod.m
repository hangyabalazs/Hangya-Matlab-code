function R = czpfs_mod(rm1,rm2)
%CZPFS   Place Field Similarity.
%   R = CZPFS_MOD(RM1,RM2) calculates linear correlation between rate maps 
%   RM1 and RM2. It discards areas where neither of the cells fire.
%
%   See also CZPLACEANALYSIS and CZPFS.

% Find non-NaNs
r1 = rm1(~isnan(rm2)&~isnan(rm1));
r2 = rm2(~isnan(rm2)&~isnan(rm1));

% Linear correlation
rr1 = r1(r1>b_max_nonnan(r1)*0.2|r2>b_max_nonnan(r2)*0.2);
rr2 = r2(r1>b_max_nonnan(r1)*0.2|r2>b_max_nonnan(r2)*0.2);
corrmtx = corrcoef(rr1,rr2);
R = corrmtx(2);