function [I_X_Y H_X H_Y I_X_Yc H_Xc H_Yc] = mutual_information(X,Y)
%MUTUAL_INFORMATION   Mutual information for discrete variables.
%   [I HX HY] = MUTUAL_INFORMATION(X,Y) calculates entropies (HX, HY) and
%   mutual information (I) for the matching samples of two discrete
%   distributions (X, Y).
%
%   [I HX HY IC HXC HYC] = MUTUAL_INFORMATION(X,Y) also returns
%   bias-corrected versions of the output variables (Panzeri-Treves bias
%   correction).
%
%   Reference:
%   Panzeri S, Senatore R, Montemurro MA, Petersen RS (2007) Correcting for
%   the sampling bias problem in spike train information measures. Journal
%   of neurophysiology 98:1064-1072
%
%   See also IMISHIFT and ENTRYRUN3C.

% Input argument check
error(nargchk(2,2,nargin))
if ~isequal(length(X),length(Y))
    error('mutual_information:inputargMismatch','Input arguments should be matching samples.')
end

% Convert the distributions to indices (integers ranging from 1)
[sx inxx] = sort(X);
[~, invx] = sort(inxx);
dx = diff(sx) > 0;
cdx = cumsum([0; dx]);
X = cdx(invx) + 1;

[sy inyy] = sort(Y);
[~, invy] = sort(inyy);
dy = diff(sy) > 0;
cdy = cumsum([0; dy]);
Y = cdy(invy) + 1;

% Joint histogram
nb1 = length(unique(X));
nb2 = length(unique(Y));
P_X_Y = accumarray([X Y],1,[nb1 nb2]);

% Marginal histograms
n = sum(P_X_Y(:));
P_X = sum(P_X_Y,2);      % X histogram
P_Y = sum(P_X_Y,1);      % Y histogram
P_X = P_X / n;       % X distribution
P_Y = P_Y / n;       % Y distribution
P_X_Y = P_X_Y / n;   % joint distribution

P_X = P_X(P_X~=0);     % restrict to the support
P_Y = P_Y(P_Y~=0);
P_X_Y = P_X_Y(P_X_Y~=0);

% Entropy
H_X = -sum(P_X.*log2(P_X));   % X ENTROPY
H_Y = -sum(P_Y.*log2(P_Y));   % Y ENTROPY
H_X_Y = -sum(P_X_Y.*log2(P_X_Y));     % COMMON ENTROPY

% Mutual information
I_X_Y = H_X + H_Y - H_X_Y;    % mutual information

% Bias correction
Nt = length(Y);   % total number of trials
R_bar = nb2;      % number of possible responses
Rs_bar = zeros(1,size(P_X_Y,1));
for k = 1:size(P_X_Y,1)
    Rs_bar(k) = bayescount(Nt,P_X_Y(k,:)/n);
end
Bias_I = (1 / (2 * Nt * log(2))) * (sum(Rs_bar-1) - (R_bar - 1));
Bias_HR = ((-1) / (2 * Nt)) * (R_bar - 1);
I_X_Yc = I_X_Y - Bias_I;
H_Xc = H_X - Bias_HR;
H_Yc = H_Y - Bias_HR;