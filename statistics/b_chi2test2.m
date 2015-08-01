function [h,p,stat] = b_chi2test2(n_d1,n_d2,alpha)
% for binned data

if nargin < 3
    alpha = 0.05;
end

ld1 = sum(n_d1);
ld2 = sum(n_d2);

% Eliminate bins that have smaller size than 4
while min(n_d1) < 4
    fnd = find(n_d1<4);
    fnd = fnd(1);
    if fnd > 1
        n_d1(fnd-1) = n_d1(fnd-1) + n_d1(fnd);
        n_d2(fnd-1) = n_d2(fnd-1) + n_d2(fnd);
    else
        n_d1(fnd+1) = n_d1(fnd+1) + n_d1(fnd);
        n_d2(fnd+1) = n_d2(fnd+1) + n_d2(fnd);
    end
    n_d1(fnd) = [];
    n_d2(fnd) = [];
end
while min(n_d2) < 4
    fnd = find(n_d2<4);
    fnd = fnd(1);
    if fnd > 1
        n_d1(fnd-1) = n_d1(fnd-1) + n_d1(fnd);
        n_d2(fnd-1) = n_d2(fnd-1) + n_d2(fnd);
    else
        n_d1(fnd+1) = n_d1(fnd+1) + n_d1(fnd);
        n_d2(fnd+1) = n_d2(fnd+1) + n_d2(fnd);
    end
    n_d1(fnd) = [];
    n_d2(fnd) = [];
end
if length(n_d1) < 2 | length(n_d2) < 2
    error('Not enough bins while discretization.')
end

stat = ld1 * ld2 * sum((n_d1/ld1-n_d2/ld2).^2./(n_d1+n_d2));
q = chi2cdf(stat,length(n_d1)-1);
if q > (1 - alpha)
    h = 1;
else
    h = 0;
end
p = 1 - q;