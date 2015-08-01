function [p, h, stats] = ranksum2(x,y)
% x >? y

alpha = 0.05;

% Get the samples and their sizes, find the larger sample
x = x(:);
y = y(:);
nx = numel(x);
ny = numel(y);
if nx <= ny
   smsample = x;
   lgsample = y;
   ns = nx;
   sg = 1;      % smsample >? lgsample
else
   smsample = y;
   lgsample = x;
   ns = ny;
   sg = 0;      % lgsample >? smsample
end

% Now deal with the method
if ns<10 && (nx+ny)<20
    method = 'exact';
else
    method = 'approximate';
end

% Compute the rank sum statistic based on the smaller sample
[ranks, tieadj] = tiedrank([smsample; lgsample]);
xrank = ranks(1:ns);
w = sum(xrank);

wmean = ns*(nx + ny + 1)/2;
if ~isequal(method,'approximate')    % use the sampling distribution of W
   if nx+ny<=10 || isequal(method,'oldexact')
      % For small samples, enumerate all possibilities, find smaller tail prob
      allpos = nchoosek(ranks,ns);
      sumranks = sum(allpos,2);
      np = length(sumranks);
      plo = sum(sumranks<=w)/np;
      phi = sum(sumranks>=w)/np;
      disp('This part needs further checking!!!')
      if sg
          p = phi;
      else
          p = plo;
      end
   else
      % Use a network algorithm for larger samples
      p = exactprob(smsample,lgsample,w);
   end
   
   p = min(p,1);           % 1-sided, p>1 means the middle is double-counted
   
else                          % use the normal approximation
   tiescor = 2 * tieadj / ((nx+ny) * (nx+ny-1));
   wvar  = nx*ny*((nx + ny + 1) - tiescor)/12;
   wc = w - wmean;
   z = (wc - 0.5 * sign(wc))/sqrt(wvar);
   p = normcdf(z,0,1);
   p = min(p,1-p);
   if (nargout > 2)
      stats.zval = z;
   end
end

if nargout > 1,
   h = (p<=alpha);
   if (nargout > 2)
      stats.ranksum = w;
   end
end

% --------------------------------------------------------------
function p = exactprob(x,y,w)
%EXACTPROB Exact P-values for Wilcoxon Mann Whitney nonparametric test
%   P=EXACTPROB(X,Y,W) computes the p-value P for the test statistic W
%   in a Wilcoxon-Mann-Whitney nonparametric test of the hypothesis that
%   X and Y come from distributions with equal medians.

% Create a contingency table with element (i,j) indicating how many
% times u(j) appears in sample i
u = unique([x(:); y(:)]);
t = zeros(2,length(u));
t(1,:) = histc(x,u)';
t(2,:) = histc(y,u)';

% Compute weights for wmw test
colsum = sum(t,1);
tmp = cumsum(colsum);
wts = [0 tmp(1:end-1)] + .5*(1+diff([0 tmp]));

% Compute p-value using network algorithm for contingency tables
p = statctexact(t,wts,w);
