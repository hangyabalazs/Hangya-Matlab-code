function vdisc = disc(unit,thres)
%DISC   Unit discriminator.
%   VD = DISC(UNIT,THR) discriminates UNIT using THR as threshold and
%   returns discriminated unit in VD.

% Input argument check
error(nargchk(2,2,nargin))
if size(unit,2) == 1
    unit = unit';
end

% Discriminating
disc0 = find(unit>=thres); 
discl = length(disc0);
disc0 = [disc0; unit(disc0(1:discl))];
disc0(1,discl+1) = length(unit) + 2;
dif = diff(disc0(1,:));
difn1 = find(dif>1);
difn1(2:length(difn1)+1) = difn1;
difn1(1) = 0;
vdisc = zeros(1,length(difn1)-1);
for j = 1:length(difn1) - 1
    [maxe,maxh] = max(disc0(2,difn1(j)+1:difn1(j+1)));
    vdisc(j) = disc0(1,difn1(j)+maxh);
end