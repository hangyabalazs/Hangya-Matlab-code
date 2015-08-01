function b_wjpth(vdtI,lentI,vdcI,lencI)

% Create 0-1 series
pstI = zeros(1,lentI);
pstI(vdtI) = 1;
pscI = zeros(1,lencI);
pscI(vdcI) = 1;

% Compress
csT = compress(pstI);
csC = compress(pscI);

% JPTH matrix
J = csT' * csC;
figure
imagesc(J);
colormap(gray)

% -------------------------------------------------------------------------
function Y = compress(X)

c = 3;
len = ceil(length(X)/c);
Y = zeros(1,len)
for i = 0:len-1
    ind1 =  c * i + 1;
    ind2 = c * (i + 1);
    ind2 = min(ind2,length(X));
    Y(i+1) = sum(X(ind1:ind2))
end