function Y = pdist2(X,str)
%PDIST2 Pairwise distance between observations.
%   Y = PDIST2(X) returns a vector which contains all the distances
%   between each pair of observations in X. X is a M by N matrix,
%   treated as M observations of N variables. Since there are 
%   M*(M-1)/2 pairs of observations in X, the size of Y is M*(M-1)/2
%   by 1. PDIST2 uses euclidean metric.
%
%   The output Y is arranged in the order of ((1,2),(1,3),..., (1,M),
%   (2,3),...(2,M),.....(M-1,M)).  i.e. the upper right triangle of
%   the M by M square matrix. To get the distance between observation
%   i and observation j, either use the formula Y((i-1)*(M-i/2)+j-i)
%   or use the helper function Z = SQUAREFORM(Y), which will return a
%   M by M symmetric square matrix, with Z(i,j) equaling the distance
%   between observation i and observation j.
%
%   Y = PDIST2(X,STR) accepts a string input taking the value 'linear' or
%   'circular' that determines the method of distance calculation. In case
%   of circular data, angles should be given in radians.
%
%   See also SQUAREFORM, LINKAGE

if nargin == 1
    str = 'linear';
end

[m, n] = size(X);

if m < 2
   error('The first argument has to be a numerical matrix with at least two rows');
end

switch str
    case 'linear'
        p = (m-1):-1:2;
        I = zeros(m*(m-1)/2,1);
        I(cumsum([1 p])) = 1;
        I = cumsum(I);
        J = ones(m*(m-1)/2,1);
        J(cumsum(p)+1) = 2-p;
        J(1)=2;
        J = cumsum(J);

        Y = zeros(1,m*(m-1)/2);
        for i = 1:m*(m-1)/2,
            K = (X(I(i),:)-X(J(i),:));
            K = sum(K.^2);
            Y(i) = sqrt(K);
        end;
    case 'circular'
        disp('Angles should be in radians.')
        Y = zeros(1,m*(m-1)/2);
        next = 1;
        for k1 = 1:m
            for k2 = k1+1:m
                Y(next) = abs(circdiff(X(k1),X(k2),'rad'));
                next = next + 1;
            end
        end
end