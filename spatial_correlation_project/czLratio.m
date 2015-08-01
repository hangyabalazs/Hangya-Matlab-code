%%

n = length(CellNumbers);
inx = find(CellNumbers==6);
cinx = find(CellNumbers~=6);

X1 = [scores_DP1(:,1)'; scores_DP1(:,2)'; scores_DP2(:,1)'; scores_DP2(:,2)'; ...
    scores_DP3(:,1)'; scores_DP3(:,2)'; scores_DP4(:,1)'; scores_DP4(:,2)'];
X = Params;
% X1 = [scores_DP1(:,1)'; scores_DP1(:,2)'; scores_DP1(:,3)'; ...
%     scores_DP2(:,1)'; scores_DP2(:,2)'; scores_DP2(:,3)'; ...
%     scores_DP3(:,1)'; scores_DP3(:,2)'; scores_DP3(:,3)'; ...
%     scores_DP4(:,1)'; scores_DP4(:,2)'; scores_DP4(:,3)'];
% % X = [X1; X];
% X1 = [scores_DP1'; scores_DP1'; scores_DP1'; ...
%     scores_DP2'; scores_DP2'; scores_DP2'; ...
%     scores_DP3'; scores_DP3'; scores_DP3'; ...
%     scores_DP4'; scores_DP4'; scores_DP4'];
% % X = [X1; X];
% z = 4;
% X1 = [scores_DP1(:,1:4)'; scores_DP2(:,1:4)'; scores_DP3(:,1:4)'; ...
%     scores_DP4(:,1:4)'];
% X = [X1; X; E1; E2; E3; E4; A1; A2; A3; A4];
% X = [X; X1]; 

XC = X(:,inx);
muC = mean(XC,2);
SC = cov(XC');
V = X - repmat(muC,1,n);
SCinv = inv(SC);

% D = zeros(1,n);
% for k = 1:n
%     D(k) = V(:,k)' * SC * V(:,k);
% end
% 
% DC = D(cinx);
% sDC = sort(DC);
% DM = sDC(length(inx));

%%

D2 = mahal(X',XC');
DC2 = D2(cinx);
sDC2 = sort(DC2);
DM2 = sDC2(length(inx))

%%

df = size(X,1);
LC = sum(1-chi2cdf(DC2,df));
nc = size(XC,2);
LrC = LC / nc