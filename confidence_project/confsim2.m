function confsim2
%CONFSIM2   Simlulate confidence based on the normative model.
%   CONFSIM2 runs a normative model for confidence and tests the pattern of
%   confidence as a function of difficulty. The external variable is
%   represented by two numbers, l (left) and r (right). The difference
%   between the numbers is randomly drawn from a discrete uniform
%   distribution of a set of numbers (representing different
%   'difficulties'). The perception model is a Gaussian distribution around
%   these numbers with fixed variance. The (deterministic) decision model
%   is a choice based on the sign of the perceived difference between the
%   numbers. 1e6 trials are simulated and confidence is based on the
%   distribution of all trials (there is no need to pre-generate a large
%   distribution and then simulate a lower number of trials, as in CONFSIM,
%   because the it is sufficient to calculate average confidence only,
%   eliminating the need for a trial-loop). Note that confidence is based
%   on the distribution of all trials, without knowledge of difficulty. The
%   results are different if confidence is based on separate distributions
%   for different difficulties. This can be simulated by changing the
%   confidence calculation - see comments in the code. Confidence is
%   calculated and plotted as a function of difficulty, separately for
%   correct and incorrect choices.
%
%   See also CONFSIM and CONFSIM3.

%   Balazs Hangya, Cold Spring Harbor Laboratory
%   1 Bungtown Road, Cold Spring Harbor
%   balazs.cshl@gmail.com
%   24-Aug-2013

%   Edit log: BH 8/24/13

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Task variables                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set default stream for pseudo-random generator
RandStream.setDefaultStream ... 
     (RandStream('mt19937ar','seed',sum(100*clock)));

% Distributions for the external variable
l_sd = 0.6;  % SD for perceived left click rate distribution (e.g. 0.6)
r_sd = 0.6;  % SD for perceived right click rate distribution 
mns = 0.01:0.25:2.51;  % different means represent different difficulities

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulate trials                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Trials
NumTrials = 1000000;  % number of trials

% Draw click rates; uniform prior
rd = randi([0 1],1,NumTrials);
rmn = randi([1 length(mns)],1,NumTrials);
l(rd==1) = 1;  % left click rate; 1 or 1+mn
l(rd==0) = 1 + mns(rmn(rd==0));
r(rd==1) = 1 + mns(rmn(rd==1));  % right click rate; 1+mn or 1
r(rd==0) = 1;
d = abs(l-r);   % difficulty

% Perceived evidence (percept) per trial
l_hat = Ldist(NumTrials,l,l_sd);  % left click rate
r_hat = Rdist(NumTrials,r,r_sd);  % right click rate
d_hat = l_hat - r_hat;  % perceived evidence ('percept')
mn = min(d_hat);
mx = max(d_hat);
bno = 100;  % bin number: d_hat has to be discretized
d_hat_inx = round((d_hat-mn)/(mx-mn)*bno);  % bin index for discretized d_hat

% Choice: sign of percept
choice = double(d_hat>0);   % choice: 1=L, 0=R

% Outcome
outcome = (l_hat > r_hat & l > r) | (l_hat < r_hat & l < r);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Confidence                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Confidence for error and correct trials
ud = unique(d_hat_inx);   % different discretized evidence/percept values
NumEVs = length(ud);   % number of diff evidence values
NumDifs = length(mns);
[confs pk0s pk1s] = deal(cell(1,NumDifs));
[pt0 pt1] = deal(zeros(1,NumDifs));   % confidence for error/correct trials for a given difficulty
for ks = 1:NumDifs   % loop through different difficluties
    S = mns(ks);  % current difficulty
    Sinx = abs(d-S) < 1e-10;   % corresponding indices
    [confs{ks} pk0s{ks} pk1s{ks}] = deal(nan(1,NumEVs));
    for k = 1:NumEVs  % loop through the evidence (d_hat) value distribution
        ct = sum(outcome==1&d_hat_inx==ud(k)) / sum(d_hat_inx==ud(k));   % confidence, calculated as outcome
%         ct = sum(outcome==1&d_hat_inx==ud(k)&Sinx) / sum(d_hat_inx==ud(k)&Sinx); % confidence, calculated as outcome, separately for different difficulties
        confs{ks}(k) = ct;
        pk0 = sum(d_hat_inx==ud(k)&outcome==0&Sinx) / sum(outcome==0&Sinx);   % probability of a given evidence on the conditional field of errors
        pk1 = sum(d_hat_inx==ud(k)&outcome==1&Sinx) / sum(outcome==1&Sinx);   % probability of a given evidence on the conditional field of correct
        pk0s{ks}(k) = pk0;
        pk1s{ks}(k) = pk1;
        pt0(ks) = pt0(ks) + nan2zero(ct*pk0);  % average confidence for errors
        pt1(ks) = pt1(ks) + nan2zero(ct*pk1);  % average confidence for correct trials
    end
end

% Plot
figure
plot(mns,pt1,'Color',[0 0.8 0],'LineWidth',4)
hold on
plot(mns,pt0,'Color',[0.8 0 0],'LineWidth',4)

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