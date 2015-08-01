function confsim2_slowfit
%CONFSIM2_SLOW   Simlulate confidence based on the normative model.
%   CONFSIM2_SLOW runs a normative model for confidence and tests the
%   psychometric performance for low and high confidence. The external
%   variable is represented by two numbers, l (left) and r (right). The
%   differentce between the numbers is randomly drawn from a discrete
%   uniform distribution of a set of numbers (representing different
%   'difficulties'). The perception model is a Gaussian distribution around
%   these numbers with fixed variance. The (deterministic) decision model
%   is a choice based on the sign of the perceived difference between the
%   numbers. Confidence is based on a pre-generated distribution of 1e6
%   trials. Simulation then runs for 10000 trials (CONFSIM2 is faster, but
%   does not generate trial-by-trial confidence values). Confidence is
%   calculated and plotted as a function of difficulty, separately for
%   correct and incorrect choices.
%
%   See also CONFSIM2.

%   Balazs Hangya, Cold Spring Harbor Laboratory
%   1 Bungtown Road, Cold Spring Harbor
%   balazs.cshl@gmail.com
%   9-Dec-2013

%   Edit log: BH 12/9/13

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Task variables                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set default stream for pseudo-random generator
RandStream.setDefaultStream ... 
     (RandStream('mt19937ar','seed',sum(100*clock)));

% Load TEmaster
load('C:\Dropbox\KepecsLab\_Josh\Data\Databases\HumanExplicitDB.mat')
inx = TEmaster.EZindex.Subject(1).AllTrials;
l = TEmaster.Precomputed.nLeftClicksExperienced_Confidence(inx);   % left clicks
r = TEmaster.Precomputed.nRightClicksExperienced_Confidence(inx);   % right clicks
d = (l-r) ./ (l+r);   % difficulty (abs of d)
mnd = min(abs(d));
mxd = max(abs(d));
bnod = 6;  % bin number: d has to be discretized
d_inx = round((abs(d)-mnd)/(mxd-mnd)*bnod);  % bin index for discretized difficulty (abs of d)
mns = unique(d_inx);  % levels of difficulty
d_sd = 0.09;

% Perceived evidence (percept) per trial
NumTrials = length(inx);   % number of trials
d_hat = Ddist(NumTrials,d,d_sd);  % perceived evidence ('percept')
mn = min(d_hat);
mx = max(d_hat);
bno = 100;  % bin number: d_hat has to be discretized
d_hat_inx = round((d_hat-mn)/(mx-mn)*bno);  % bin index for discretized d_hat

% Choice: sign of percept
choice = double(d_hat>0);   % choice: 1=L, 0=R

% Outcome
outcome = (d_hat > 0 & l > r) | (d_hat < 0 & l < r);
oinx = find(l==r);
outcome(oinx) = randi([0 1],size(oinx));   % random outcome for equal rates

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulate confidence                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Trials
conf = nan(1,NumTrials);
for iT = 1:NumTrials   % loop through trials
    
    % Confidence per trial
    if choice(iT)   % left side choice
        P_H_L_d_hat = sum(l>r&d_hat_inx==d_hat_inx(iT)) / NumTrials;   % P(H_L,d_hat)
        P_d_hat = sum(d_hat_inx==d_hat_inx(iT)) / NumTrials;   % P(d_hat)
        conf(iT) = P_H_L_d_hat / P_d_hat;   % confidence per trial
    else   % right side choice
        P_H_R_d_hat = sum(r>l&d_hat_inx==d_hat_inx(iT)) / NumTrials;   % P(H_R,d_hat)
        P_d_hat = sum(d_hat_inx==d_hat_inx(iT)) / NumTrials;   % P(d_hat)
        conf(iT) = P_H_R_d_hat / P_d_hat;   % confidence per trial
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Confidence vs. percept                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Confidence for error and correct trials
NumDifs = length(mns);
[pt0 pt1] = deal(zeros(1,NumDifs));   % confidence for error/correct trials for a given difficulty
for ks = 1:NumDifs   % loop through different difficluties
    S = mns(ks);  % current difficulty
    Sinx = abs(abs(d_inx)-S) < 1e-10;   % corresponding indices
    pt0(ks) = mean(conf(Sinx&outcome==0));   % confidence for error trials with current difficulty level
    pt1(ks) = mean(conf(Sinx&outcome==1));   % confidence for correct trials with current difficulty level
end

% Plot
figure
plot(mns,pt1,'Color',[0 0.8 0],'LineWidth',4)
hold on
plot(mns,pt0,'Color',[0.8 0 0],'LineWidth',4)

keyboard

% -------------------------------------------------------------------------
function D = Ddist(N,M,SD)

D = SD * randn(1,N) + M;