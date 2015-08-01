function [h,p,stat] = b_chi2test(d1,d2,alpha)

if nargin < 3
    alpha = 0.05;
end

ld1 = length(d1);
ld2 = length(d2);

% Discreatization
% hb = fix(exp(0.626+0.4*log(min(length(d1),length(d2)))));     % number of bins
hb = round(min(ld1,ld2)/4);
mx = max([d1 d2]);
mn = min([d1 d2]);
pmx = chi2cdf(mx,hb-1);
pmn = chi2cdf(mn,hb-1);
inve = linspace(pmn,pmx,hb);
edges = chi2inv(inve,hb-1);     % Más cdf inverze kéne!!!
% min = min(min([d1 d2]));
% max = max(max([d1 d2]));
% binwidth = (max - min) / hb;
% edges = [min:binwidth:max];
n_d1 = histc(d1,edges);
n_d1(end-1) = n_d1(end-1) + n_d1(end);
n_d1 = n_d1(1:end-1);
n_d2 = histc(d2,edges);
n_d2(end-1) = n_d2(end-1) + n_d2(end);
n_d2 = n_d2(1:end-1);

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