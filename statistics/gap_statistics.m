function [khat Gap s_k C] = gap_statistics(X,K,method)
%GAP_STATISTICS   Optimal number of clusters.
%   KHAT = GAP_STATISTICS(X,K,METHOD) performs hierarchical clustering on X
%   using the algorithm passed in METHOD (optional; default: 'single', i.e
%   shortest distance) and squared Euclidean distance. If X is
%   3-dimensional, distance is calculated for each 2D matrices and
%   averaged. It calcluates the gap statistics (see reference below) from
%   number of clusters ranging from 1 to K (optional; default: 15).
%   Bootstrap for within-cluster sum of squares is calculated using uniform
%   reference distributions over a box aligned to the principal components
%   of the data ('Gap/pc' method), with 20 Monte Carlo replicates per
%   number of clusters. Optimal number of clusters is returned in KHAT. The
%   gap statistics is plotted against the number of clusters.
%
%   [KHAT GAP S_K] = GAP_STATISTICS(X,K,METHOD) returns the gap statistics
%   (GAP) and the standard deviation of the bootstrap sample (S_K), too.
%
%   [KHAT GAP S_K C] = GAP_STATISTICS(X,K,METHOD) returns the cluster
%   assignment matrix C (see CLUSTER).
%
%   Reference:
%   Tibshiriani R, Walther G, Hastie T (2001) Estimating the number of
%   clusters in a data set via the gap statistic. J R Statist Soc B
%   63:411-423
%
%   See also SOMCCGCLUST2.

% Input arguments check
error(nargchk(1,3,nargin));
if nargin < 3
    method = 'single';     % set default method to 'single' - shortest ditance
end
if nargin < 2
    K = 15;    % default max. cluster number
end

% Cluster the data; calculate within-cluster dispersion
[W C] = gapcluster(X,K,method);

% Reference distribution
B = 20;   % bootstrap size (Monte Carlo replicates)
Wstar = nan(B,K);
for b = 1:B
    dims = length(size(X));
    if dims > 2   % 3D input matrix
        dds = size(X,3);
        Z = nan(size(X));
        for k = 1:dds
            Xk = squeeze(X(:,:,k));
            Z(:,:,k) = refdistr(Xk);    % uniform reference distribution over a box aligned with the principal components of the data
        end
    else
        Z = refdistr(X);    % uniform reference distribution over a box aligned with the principal components of the data
    end
    
    % Cluster the reference distribution
    Wstar(b,:) = gapcluster(Z,K,method);   % bootstrap within-cluster dispersion
end

% Gap statistics
Gap = (1/B) * sum(log(Wstar)-repmat(log(W),B,1));

% Standard deviation
% lbar = (1/B) * sum(log(Wstar));   % equivalent calculation of the SD of log(Wstar)
% sd_k = sqrt((1/B)*sum((log(Wstar)-repmat(lbar,B,1)).^2));
sd_k = std(log(Wstar),1);
s_k = sd_k * sqrt(1+1/B);   % correct for simulation error

% Number of clusters
lg = Gap(1:end-1) >= Gap(2:end) - s_k(2:end);
khat = find(lg,1,'first');      % smallest k such that Gap(k)>=Gap(k+1)-s(k+1)

% Visualize
figure
plot(1:K,Gap,'k')
hold on
errorbar(1:K,Gap,s_k,'k+')

% -------------------------------------------------------------------------
function [W C] = gapcluster(X,K,method)

% Pairwise distances
dims = length(size(X));
if dims > 2   % 3D input matrix
    dds = size(X,3);     % calculate distance for all input variables (3rd dim)
    pds = nan(dds,nchoosek(size(X,1),2));
    for k = 1:dds
        pds(k,:) = pdist(squeeze(X(:,:,k))) .^ 2;   % squared Euclidean distance
    end
    dst = mean(pds);   % average over input variables
else
    dst = pdist(X) .^ 2;    % one variable (traditional); squared Euclidean distance
end
links = linkage(dst,method);

% Iterative clustering; calculate within-cluster dispersion
qdst = squareform(dst);     % square format distance matrix
W = nan(1,K);
C = cluster(links,'maxclust',1:K);   % iterative clustering over different number of clusters
for k = 1:K   % within-cluster dispersion over different number of clusters
    c = C(:,k);
    D = nan(1,k);
    Dnorm = nan(1,k);
    for r = 1:k    % loop through all clusters to determine within cluster distances
        inxr = find(c==r);  % elements of the rth cluster
        n_r = length(inxr); % number of elements in the cluster
        d_ij = qdst(inxr,inxr);     % within cluster distances
        D(r) = sum(d_ij(:));    % sum of within cluster distances
        Dnorm(r) = D(r) / (2 * n_r);    % equivalent to within-cluster sum of squares for squared Euclidean distance
    end
    W(k) = sum(Dnorm);   % within-cluster dispersion
end

% -------------------------------------------------------------------------
function Z = refdistr(X)

% Reference distribution
mX = mean(X);    % mean of columns
n = size(X,1);   % number of points
Xnorm = X - repmat(mX,n,1);  % transform to 0 mean
[U,S,V] = svd(Xnorm);   %#ok<ASGLU> % singular value decomposition
Xprime = Xnorm * V';
Zprime = uniformref(Xprime);    % uniform reference distribution over a box aligned with the principal components of the data
Z = Zprime * V;     % transform back to data space
Z = Z + repmat(mX,n,1);     % transform back to original mean

% -------------------------------------------------------------------------
function Z = uniformref(X)

% Input range
mnX = min(X);
mxX = max(X);

% Uniform distribution over the range
n = size(X,1);
Z = rand(size(X)) .* repmat(mxX-mnX,n,1) + repmat(mnX,n,1);