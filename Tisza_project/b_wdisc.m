function [vdtI vdcI] = b_wdisc(tI,cI)
%WDISC   Discriminate Tisza data.
%
%   See also WREAD and DISC.

thres1 = 50;
vdtI = discrim(tI,thres1);
thres2 = 50;
vdcI = discrim(cI,thres2);

% -------------------------------------------------------------------------
function vdisc = discrim(unit,kuszob)

disc = find(unit>=kuszob);
discl = length(disc);
disc = [disc; unit(disc(1:discl))];
disc(1,discl+1) = length(unit) + 2;
dif = diff(disc(1,:));
difn1 = find(dif>1);
difn1(2:length(difn1)+1) = difn1;
difn1(1) = 0;
vdisc = zeros(1,length(difn1)-1);
for j = 1:length(difn1) - 1
    [maxe,maxh] = max(disc(2,difn1(j)+1:difn1(j+1)));
    vdisc(j) = disc(1,difn1(j)+maxh);
end