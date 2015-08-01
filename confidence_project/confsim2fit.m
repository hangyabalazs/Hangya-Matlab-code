function confsim2fit
%CONFSIM2FIT   Fit the normative model.
%   CONFSIM2FIT uses a psychometric fit on human performance data to model
%   perception in the normative confidence model (see CONFSIM and CONFSIM2
%   for model details). It loads the human confidence dataset and runs the
%   model on the presented evidence streams with the Gaussain noise model
%   of perception using the fitted psychometric parameter as standard
%   deviation. It plot confidence for correct and error trials as a
%   function of difficulty. Difficulty is defined as absolute difference
%   over sum of the presented eveidence in each trial.
%
%   See also CONFSIM and CONFSIM2.

%   Balazs Hangya, Cold Spring Harbor Laboratory
%   1 Bungtown Road, Cold Spring Harbor
%   balazs.cshl@gmail.com
%   19-Sept-2013

%   Edit log: BH 9/19/13

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Task variables                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set default stream for pseudo-random generator
RandStream.setDefaultStream ... 
     (RandStream('mt19937ar','seed',sum(100*clock)));

% Load TEmaster
load('C:\Dropbox\KepecsLab\_Josh\Data\Databases\HumanExplicitDB.mat')
inx = TEmaster.EZindex.Subject(1).AllTrials;
l = TEmaster.Precomputed.nLeftClicksExperienced_Confidence(inx);
r = TEmaster.Precomputed.nRightClicksExperienced_Confidence(inx);
d = (l-r)./(l+r);
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
% Confidence                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Confidence for error and correct trials
ud = unique(d_hat_inx);   % different discretized evidence/percept values
NumEVs = length(ud);   % number of diff evidence values
NumDifs = length(mns);  % number of tril difficulties
[confs pk0s pk1s] = deal(cell(1,NumDifs));
[pt0 pt1] = deal(zeros(1,NumDifs));   % confidence for error/correct trials for a given difficulty
for ks = 1:NumDifs   % loop through different difficluties
    S = mns(ks);  % current difficulty
    Sinx = abs(abs(d_inx)-S) < 1e-10;   % corresponding indices
    [confs{ks} pk0s{ks} pk1s{ks}] = deal(nan(1,NumEVs));
    for k = 1:NumEVs  % loop through the evidence (d_hat) value distribution
        ct = sum(outcome==1&d_hat_inx==ud(k)) / sum(d_hat_inx==ud(k));   % confidence, calculated as outcome
%         ct = sum(outcome==1&d_hat_inx==ud(k)&Sinx) / sum(d_hat_inx==ud(k)&Sinx); % confidence, alculated as outcome, separately for different difficulties
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
function D = Ddist(N,M,SD)

D = SD * randn(1,N) + M;