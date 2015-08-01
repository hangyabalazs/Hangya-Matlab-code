function confsim2_slow
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

% Generate series (distribution) of percepts
N = 1000000;
L_SD = 0.6;  % SD for perceived left click rate distribution (e.g. 0.6)
R_SD = 0.6;  % SD for perceived right click rate distribution 
MNs = 0.01:0.25:2.51;  % different means represent different difficulities
Rd = randi([0 1],1,N);   % left/right
Rmn = randi([1 length(MNs)],1,N);   % difficulty levels
L(Rd==1) = 1;  % left click rate; 1 or 1+MN
L(Rd==0) = 1 + MNs(Rmn(Rd==0));
R(Rd==1) = 1 + MNs(Rmn(Rd==1));  % right click rate; 1+MN or 1
R(Rd==0) = 1;
L_hat = Ldist(N,L,L_SD);  % percieved left click rate
R_hat = Rdist(N,R,R_SD);  % percieved right click rate
D_hat = L_hat - R_hat;  % perceived evidence ('percept')
Perf = sum((D_hat>0&L>R)|(D_hat<0&L<R)) / N;   % average performance

% Discretize percept
mn = min(D_hat);
mx = max(D_hat);
bno = 100;  % bin number: d_hat has to be discretized
D_hat_disc = round((D_hat-mn)/(mx-mn)*bno);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulate trials                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Trials
NumTrials = 10000;  % number of trials

% Draw click rates; uniform prior
rd = randi([0 1],1,NumTrials);   % left/right
rmn = randi([1 length(MNs)],1,NumTrials);   % difficulty levels
l(rd==1) = 1;  % left click rate; 1 or 1+mn
l(rd==0) = 1 + MNs(rmn(rd==0));
r(rd==1) = 1 + MNs(rmn(rd==1));  % right click rate; 1+mn or 1
r(rd==0) = 1;
d = abs(l-r);   % difficulty

% Perceived evidence (percept) per trial
l_hat = Ldist(NumTrials,l,L_SD);  % left click rate
r_hat = Rdist(NumTrials,r,R_SD);  % right click rate
d_hat = l_hat - r_hat;  % perceived evidence ('percept')
d_hat_inx = round((d_hat-mn)/(mx-mn)*bno);  % bin index for discretized d_hat

% Choice: sign of percept
choice = double(d_hat>0);   % choice: 1=L, 0=R

% Outcome
outcome = (l_hat > r_hat & l > r) | (l_hat < r_hat & l < r);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Confidence                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Trials
conf = nan(1,NumTrials);
for iT = 1:NumTrials   % loop through trials
    
    % Confidence per trial
    if choice(iT)   % left side choice
        P_H_L_d_hat = sum(L>R&D_hat_disc==d_hat_inx(iT)) / N;   % P(H_L,d_hat)
        P_d_hat = sum(D_hat_disc==d_hat_inx(iT)) / N;   % P(d_hat)
        conf(iT) = P_H_L_d_hat / P_d_hat;   % confidence per trial
    else   % right side choice
        P_H_R_d_hat = sum(R>L&D_hat_disc==d_hat_inx(iT)) / N;   % P(H_R,d_hat)
        P_d_hat = sum(D_hat_disc==d_hat_inx(iT)) / N;   % P(d_hat)
        conf(iT) = P_H_R_d_hat / P_d_hat;   % confidence per trial
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Confidence vs. percept                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Confidence for error and correct trials
NumDifs = length(MNs);
[pt0 pt1] = deal(zeros(1,NumDifs));   % confidence for error/correct trials for a given difficulty
for ks = 1:NumDifs   % loop through different difficluties
    S = MNs(ks);  % current difficulty
    Sinx = abs(d-S) < 1e-10;   % corresponding indices
    pt0(ks) = mean(conf(Sinx&outcome==0));   % confidence for error trials with current difficulty level
    pt1(ks) = mean(conf(Sinx&outcome==1));   % confidence for correct trials with current difficulty level
end

% Plot
figure
plot(MNs,pt1,'Color',[0 0.8 0],'LineWidth',4)
hold on
plot(MNs,pt0,'Color',[0.8 0 0],'LineWidth',4)

keyboard

% -------------------------------------------------------------------------
function D = Ldist(N,M,SD)

D = SD * randn(1,N) + M;
% r = randi([0 1],1,N);
% D = r .* (SD .* randn(1,N) + M) + (1 - r) .* (SD / 2 * randn(1,N) + M - 2);

% -------------------------------------------------------------------------
function D = Rdist(N,M,SD)

D = SD * randn(1,N) + M;
% r = randi([0 1],1,N);
% D = r .* (SD .* randn(1,N) + M) + (1 - r) .* (SD / 2 * randn(1,N) + M - 2);