function confsim3b
%CONFSIM3B   Simlulate confidence based on the normative model.
%   CONFSIM3B runs a normative model for confidence and tests the
%   psychometric performance for low and high confidence. The external
%   variable is represented by two numbers, l (left) and r (right). The
%   differentce between the numbers is randomly drawn from a discrete
%   uniform distribution of a set of numbers (representing different
%   'difficulties'). The perception model is a Gaussian distribution around
%   these numbers with fixed variance. The (deterministic) decision model
%   is a choice based on the sign of the perceived difference between the
%   numbers. 10000 trials are simulated and confidence is based on the
%   distribution of all trials (this gives a low resolution distribution
%   for confidence calculation, resulting in less CPU time but noisier
%   curves; for slower but smoother results, see CONFSIM3). Note that
%   confidence is based on the distribution of all trials, without
%   knowledge of difficulty. Psychometric functions are plotted for low and
%   high confidence trials (median split).
%
%   See also CONFSIM2, CONFSIM3, CONFSIM3 and CONFSIM3C.

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
NumTrials = 10000;  % number of trials

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
% Performance depend on confidence                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Psychometric performance for low and high confidence trials
NumDifs = length(mns);
[Perf_lowC Perf_highC] = deal(zeros(1,NumDifs));
for ks = 1:NumDifs   % loop through different difficluties
    S = mns(ks);  % current difficulty
    Sinx = abs(d-S) < 1e-10;   % corresponding indices
    Perf_lowC(ks) = mean(outcome(Sinx&conf<median(conf)));   % low confidence trials
    Perf_highC(ks) = mean(outcome(Sinx&conf>median(conf)));   % high confidence trials
end

% Plot
figure
plot(mns,Perf_lowC,'Color',[0.2039 0.3020 0.4941],'LineWidth',4)
hold on
plot(mns,Perf_highC,'Color',[0 0.749 0.749],'LineWidth',4)

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