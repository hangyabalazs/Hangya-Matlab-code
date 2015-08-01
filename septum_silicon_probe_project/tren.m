function [TE TE_corr TE_bias H_X2FcX2] = tren(W1,W2,sr,tau)
%TREN   Transfer entropy.
%   TE = TREN(W1,W2,SR) calculates transfer entropy for time series W1 and
%   W2 sampled on SR (TE_W1->W2).
%
%   TE = TREN(W1,W2,SR,TAU) uses TAU ms as time lag between current and
%   "future" values.
%
%   [TE TE_CORR TE_BIAS] = TREN(W1,W2,SR,TAU) returns unbiased estimate of
%   transfer entropy applying Treves-Panzeri bias correction algorithm. The
%   amount of the bias of the original estimate is also returned.
%
%   [TE TE_CORR TE_BIAS H] = TREN(W1,W2,SR,TAU) returns H(X2F|X2)
%   conditonal entropy for transfer entropy normalization.
%
%   References:
%   (1) Panzeri S, Senatore R, Montemurro MA, Petersen RS (2007)
%   Correcting for the sampling bias problem in spike train information
%   measures. J Neurophysiol 98:1064-1072.
%   (2) Gourévitch B, Eggermont JJ (2007) Evaluating information transfer
%   between auditory cortical neurons. J Neurophysiol 97:2533-2543.
%   (3) Imas OA, Ropella KM, Ward BD, Wood JD, Hudetz AG (2005) Volatile
%   anesthetics disrupt frontal-posterior recurrent information transfer at
%   gamma frequencies in rat. Neurosci Lett 387:145-50.

% Input argument check
error(nargchk(3,4,nargin))
if nargin < 4
    tau = 100;   % time lag in ms
end    

% Extracting the random variables
tau2 = tau / 1000 * sr;     % time lag in data points
X1 = W1(1:end-tau2);
X2 = W2(1:end-tau2);
X2F = W2(tau2:end+1);
n = length(X1);

% Calculating joint histograms
bno = fix(exp(0.626+0.4*log(n-1)));   % bin number for histogram estimates
bno = 30;
minX1 = min(X1);    % discretization of X1
maxX1 = max(X1);
binwidth1 = (maxX1 - minX1) ./ bno;
xx1 = minX1 + binwidth1 * (0:bno);   % bin edges
xx1(length(xx1)) = maxX1;
xx1(1) = -inf;
nbin1 = length(xx1);

minX2 = min(W2);    % discretization of X2
maxX2 = max(W2);
binwidth2 = (maxX2 - minX2) ./ bno;
xx2 = minX2 + binwidth2 * (0:bno);   % bin edges
xx2(length(xx2)) = maxX2;
xx2(1) = -inf;
nbin2 = length(xx2);

h_X2F_X2_X1 = zeros(nbin2-1,nbin2-1,nbin1-1);   % P(X2F, X2, X1)
t1 = X2F(:) - minX2;
t2 = X2(:) - minX2;
t3 = X1(:) - minX1;
p1 = fix((t1-1000000*eps)/binwidth2) + 1;
p2 = fix((t2-1000000*eps)/binwidth2) + 1;
p3 = fix((t3-1000000*eps)/binwidth1) + 1;
for k = 1:n
    ind1 = min(p1(k),size(h_X2F_X2_X1,1));
    ind2 = min(p2(k),size(h_X2F_X2_X1,2));
    ind3 = min(p3(k),size(h_X2F_X2_X1,3));
    h_X2F_X2_X1(ind1,ind2,ind3) = h_X2F_X2_X1(ind1,ind2,ind3) + 1;
end
h_X2F_X2_X1 = h_X2F_X2_X1 / sum(sum(sum(h_X2F_X2_X1)));      % normalization

% h_X2_X1 = zeros(nbin2-1,nbin1-1);   % P(X2, X1)
% t1 = X2(:) - minX2;
% t2 = X1(:) - minX1;
% p1 = fix((t1-1000000*eps)/binwidth2) + 1;
% p2 = fix((t2-1000000*eps)/binwidth1) + 1;
% for k = 1:n
%     ind1 = min(p1(k),size(h_X2_X1,1));
%     ind2 = min(p2(k),size(h_X2_X1,2));
%     h_X2_X1(ind1,ind2) = h_X2_X1(ind1,ind2) + 1;
% end
% 
% h_X2F_X2 = zeros(nbin2-1,nbin2-1);   % P(X2F, X2)
% t1 = X2F(:) - minX2;
% t2 = X2(:) - minX2;
% p1 = fix((t1-1000000*eps)/binwidth2) + 1;
% p2 = fix((t2-1000000*eps)/binwidth2) + 1;
% for k = 1:n
%     ind1 = min(p1(k),size(h_X2F_X2,1));
%     ind2 = min(p2(k),size(h_X2F_X2,2));
%     h_X2F_X2(ind1,ind2) = h_X2F_X2(ind1,ind2) + 1;
% end

% Calculating marginal histograms
h_X2_X1 = squeeze(sum(h_X2F_X2_X1,1));   % P(X2, X1)
h_X2F_X2 = squeeze(sum(h_X2F_X2_X1,3));   % P(X2F, X2)
h_X2 = squeeze(sum(h_X2_X1,2));          % P(X2)

% Calculating transfer entropy
TE = 0;
for k1 = 1:nbin2-1
    for k2 = 1:nbin2-1
        for k3 = 1:nbin1-1
            if h_X2F_X2_X1(k1,k2,k3) ~= 0
                Pa = h_X2F_X2_X1(k1,k2,k3);
                Pb = h_X2(k2);
                Pc = h_X2_X1(k2,k3);
                Pd = h_X2F_X2(k1,k2);
                TE = TE + Pa * log2((Pa*Pb)/(Pc*Pd));
            end
        end
    end
end

% Alternative calculation of transfer entropy (for bias correction purposes)
H_X2F_X2 = 0;       % H(X2F, X2)
for k1 = 1:nbin2-1
    for k2 = 1:nbin2-1
        if h_X2F_X2(k1,k2) ~= 0
            Pa = h_X2F_X2(k1,k2);
            H_X2F_X2 = H_X2F_X2 - Pa * log2(Pa);
        end
    end
end
H_X2 = 0;           % H(X2)
for k2 = 1:nbin2-1
    if h_X2(k2) ~= 0
        Pa = h_X2(k2);
        H_X2 = H_X2 - Pa * log2(Pa);
    end
end
H_X2FcX2 = H_X2F_X2 - H_X2;     % H(X2F|X2)

H_X2F_X2_X1 = 0;    % H(X2F, X2, X1)
for k1 = 1:nbin2-1
    for k2 = 1:nbin2-1
        for k3 = 1:nbin1-1
            if h_X2F_X2_X1(k1,k2,k3) ~= 0
                Pa = h_X2F_X2_X1(k1,k2,k3);
                H_X2F_X2_X1 = H_X2F_X2_X1 - Pa * log2(Pa);
            end
        end
    end
end
H_X2_X1 = 0;    % H(X2,X1)
for k2 = 1:nbin2-1
    for k3 = 1:nbin1-1
        if h_X2_X1(k2,k3) ~= 0
            Pa = h_X2_X1(k2,k3);
            H_X2_X1 = H_X2_X1 - Pa * log2(Pa);
        end
    end
end
H_X2FcX2_X1 = H_X2F_X2_X1 - H_X2_X1;    % H(X2F|X1,X2)
TE2 = H_X2FcX2 - H_X2FcX2_X1;

% Bias correction
Nt = n;   % total number of trials
Rs_bar = zeros(1,size(h_X2F_X2,2));
for k2 = 1:size(h_X2F_X2,2)
    Rs_bar(k2) = bayescount(Nt,h_X2F_X2(:,k2));
end
Bias_HRS = ((-1) / (2 * Nt * log(2))) * sum(Rs_bar-1);
H_X2FcX2_corr = H_X2FcX2 - Bias_HRS;

Nt = n;   % total number of trials
Rs_bar = zeros(size(h_X2F_X2_X1,2),size(h_X2F_X2_X1,3));
for k2 = 1:size(h_X2F_X2_X1,2)
    for k3 = 1:size(h_X2F_X2_X1,3)
        Rs_bar(k2,k3) = bayescount(Nt,h_X2F_X2_X1(:,k2,k3));
    end
end
Bias_HRS = ((-1) / (2 * Nt * log(2))) * sum(sum(Rs_bar-1));
H_X2FcX2_X1_corr = H_X2FcX2_X1 - Bias_HRS;

TE_corr = H_X2FcX2_corr - H_X2FcX2_X1_corr;
TE_bias = TE - TE_corr;