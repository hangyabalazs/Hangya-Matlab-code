function confsim2_slowfit_fit_joshedit(TEmaster, SubjectID, nEvidenceBins, PerceptualNoiseFactor, ConfidenceNoiseFactor, nBoot)
%CONFSIM2_SLOW   Simlulate confidence based on the normative model.
%   CONFSIM2_SLOW runs a normative model for confidence and tests the
%   psychometric performance for low and high confidence. The external
%   variable is represented by two numbers, l (left) and r (right). The
%   differentce between the numbers is randomly drawn from a discrete
%   uniform distribution of a set of numbers (representing different
%   'difficulties'). The perception model is a Gaussian distribution around
%   these numbers with fixed variance. The (deterministic) decision model
%   is a choice based on the sign of the perceived difference between the
%   numbers. Confidence is 4based on a pre-generated distribution of 1e6
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

if matlabpool('size') == 0
    matlabpool local 7
end

% Set default stream for pseudo-random generator
RandStream.setGlobalStream ... 
     (RandStream('mt19937ar','seed',sum(100*clock)));

% load('C:\Dropbox\KepecsLab\_Josh\Data\Databases\HumanExplicitDB.mat')
SubjectTrials = [];
for x = 1:length(SubjectID)
    SubjectTrials = [SubjectTrials TEmaster.EZindex.Subject(SubjectID(x)).AllTrials]; 
end
inx = SubjectTrials;
inx = intersect(inx, find(TEmaster.ReactionTime > .05)); % Include only trials with sampled evidence
inx = intersect(inx, find(TEmaster.ChosenDirection < 3)); % Include only completed trials

l = TEmaster.Precomputed.nLeftClicksExperienced_Confidence(inx);   % left clicks
r = TEmaster.Precomputed.nRightClicksExperienced_Confidence(inx);   % right clicks
so = TEmaster.Precomputed.OutcomeFromExperiencedRates(inx); % human outcomes
d = (l-r) ./ (l+r);   % difficulty (abs of d)
mnd = min(abs(d));
mxd = max(abs(d));
bnod = nEvidenceBins;  % bin number: d has to be discretized
d_inx = round((abs(d)-mnd)/(mxd-mnd)*(bnod-1))+1;  % bin index for discretized difficulty (abs of d) - corrected so there isn't a zero index. JS
mns = unique(d_inx);  % levels of difficulty
d_sd = PerceptualNoiseFactor;

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

% Transfer function (discretized difficulty -> distribution of confidence)
[confTF outcomeTF] = deal(nan(bnod,nBoot));
parfor d_inxTF = 1:bnod
    allds = d(d_inx==d_inxTF);
%     dTF = sign(rand(1)-0.5) * d_inxTF / bnod * (mxd - mnd) + mnd;   % mid-percept from discretized percept
    dTF = allds(randi([1 length(allds)],1,1));   % draw corresponding difficulty from the distribution underlying the discretized diff. in the data
    for k = 1:nBoot
        d_hatTF = Ddist(1,dTF,d_sd);   % percept
        d_hat_inxTF = round((d_hatTF-mn)/(mx-mn)*bno);  % discretized percept
        choiceTF = double(d_hatTF>0);   % choice: 1=L, 0=R
        if choiceTF   % left side choice
            P_H_L_d_hat = sum(l>r&d_hat_inx==d_hat_inxTF) / NumTrials;   % P(H_L,d_hat)
            P_d_hat = sum(d_hat_inx==d_hat_inxTF) / NumTrials;   % P(d_hat)
            confTF(d_inxTF,k) = P_H_L_d_hat / P_d_hat;   % confidence per trial
        else   % right side choice
            P_H_R_d_hat = sum(r>l&d_hat_inx==d_hat_inxTF) / NumTrials;   % P(H_R,d_hat)
            P_d_hat = sum(d_hat_inx==d_hat_inxTF) / NumTrials;   % P(d_hat)
            confTF(d_inxTF,k) = P_H_R_d_hat / P_d_hat;   % confidence per trial
        end
        outcomeTF(d_inxTF,k) = (d_hatTF > 0 & dTF > 0) | (d_hatTF < 0 & dTF < 0);
    end
end

% Compute confidence lookup table for correct and error
CorrectConfTF = zero2nan(confTF.*outcomeTF);
ErrorConfTF = zero2nan(confTF.*(1-outcomeTF));
CorrectConfTFMean = nanmean(CorrectConfTF');
ErrorConfTFMean = nanmean(ErrorConfTF');

% Filter means with insufficient data
SufficiencyCutoff = 30; % at least 30 points contributing to mean
CorrectCounts = sum(~isnan(CorrectConfTF'));
ErrorCounts = sum(~isnan(ErrorConfTF'));
CorrectConfTFMean(CorrectCounts < SufficiencyCutoff) = NaN;
ErrorConfTFMean(ErrorCounts < SufficiencyCutoff) = NaN;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulate confidence                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generate predicted confidence for each trial
conf = nan(1,NumTrials);
for x = 1:NumTrials
    if so(x) == 1
        conf(x) = CorrectConfTFMean(d_inx(x));
    elseif so(x) == 0
        conf(x) = ErrorConfTFMean(d_inx(x));
    end
end

ModelConf_binned = Map2FivePtScale(conf);
ModelConf_binned(ModelConf_binned==0) = NaN;

% parfor iT = 1:NumTrials   % loop through trials
%     
%     % Confidence per trial
%     if choice(iT)   % left side choice
%         P_H_L_d_hat = sum(l>r&d_hat_inx==d_hat_inx(iT)) / NumTrials;   % P(H_L,d_hat)
%         P_d_hat = sum(d_hat_inx==d_hat_inx(iT)) / NumTrials;   % P(d_hat)
%         conf(iT) = P_H_L_d_hat / P_d_hat;   % confidence per trial
%     else   % right side choice
%         P_H_R_d_hat = sum(r>l&d_hat_inx==d_hat_inx(iT)) / NumTrials;   % P(H_R,d_hat)
%         P_d_hat = sum(d_hat_inx==d_hat_inx(iT)) / NumTrials;   % P(d_hat)
%         conf(iT) = P_H_R_d_hat / P_d_hat;   % confidence per trial
%     end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Confidence vs. percept                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Confidence for error and correct trials
NumDifs = length(mns);
[pt0 pt1] = deal(zeros(1,NumDifs));   % confidence for error/correct trials for a given difficulty
for ks = 1:NumDifs   % loop through different difficluties
    S = mns(ks);  % current difficulty
    Sinx = abs(abs(d_inx)-S) < 1e-10;   % corresponding indices
    pt0(ks) = mean(conf(Sinx&so==0));   % confidence for error trials with current difficulty level
    pt1(ks) = mean(conf(Sinx&so==1));   % confidence for correct trials with current difficulty level
end

% Plot
figure
plot(mns-1,pt1,'Color',[0 0.8 0],'LineWidth',4)
hold on
plot(mns-1,pt0,'Color',[0.8 0 0],'LineWidth',4)

% Compute human confidence for same evidence bins, mapped to .5-1 scale
HumanOutcomes = TEmaster.Precomputed.OutcomeFromExperiencedRates(inx);
HumanConfidenceReports = TEmaster.ExplicitConfidenceReport(inx);
[CorrectConfidence ErrorConfidence CorrectConfidenceError ErrorConfidenceError] = deal(nan(1,nEvidenceBins));
SubjectMinPoints = 10;
for x = mns
    CorrectConfidenceReports = HumanConfidenceReports(d_inx==x&HumanOutcomes==1);
    %CorrectConfidenceReports = ((CorrectConfidenceReports/5)/2)+.45; % Map to .5 --> 1 scale with even bins (1=.55 2=.65,etc)
    CorrectConfidenceReports = Rail2RailMap(CorrectConfidenceReports);
    CorrectConfidence(x) = nanmean(CorrectConfidenceReports);
    CorrectConfidenceError(x) = serr(CorrectConfidenceReports)*1.96;
    ErrorConfidenceReports = HumanConfidenceReports(d_inx==x&HumanOutcomes==0);
    %ErrorConfidenceReports = ((ErrorConfidenceReports/5)/2)+.45;
    ErrorConfidenceReports = Rail2RailMap(ErrorConfidenceReports);
    ErrorConfidence(x) = nanmean(ErrorConfidenceReports);
    ErrorConfidenceError(x) = serr(ErrorConfidenceReports)*1.96;
    if length(CorrectConfidenceReports) < SubjectMinPoints
        CorrectConfidence(x) = NaN;
    end
    if length(ErrorConfidenceReports) < SubjectMinPoints
        ErrorConfidence(x) = NaN;
    end
end
errorbar(mns-1,CorrectConfidence, CorrectConfidenceError, 'Color', [0 0.8 0], 'LineWidth', 4, 'LineStyle', 'none'); hold on;
errorbar(mns-1,ErrorConfidence, ErrorConfidenceError, 'Color', [1 0 0], 'LineWidth', 4, 'LineStyle', 'none');
k = 5;
set(gca, 'Box', 'off', 'TickDir', 'out', 'FontName', 'Arial', 'TickLength', [0.025 0.05], 'FontSize', 18);
set(gca, 'YLim', [.5 1], 'YTick', [.5 .6 .7 .8 .9 1], 'XLim', [0 length(mns)])
ylabel(gca,'Confidence (pCorrect)', 'FontSize', 18, 'FontName', 'Arial');

% Calibration plot
figure;
SubjectOutcomesBinned = cell(1,5); ModelOutcomesBinned = cell(1,5); % Cell array of original outcomes in each bin in case its needed later for fitting, etc
Mean = nan(1,5); ModelMean = nan(1,5);
Error = nan(1,5); ModelError = nan(1,5);
Outcomes = TEmaster.Precomputed.OutcomeFromExperiencedRates(SubjectTrials);
Confidence = TEmaster.ExplicitConfidenceReport(SubjectTrials);
for x = 1:5
    SubjectOutcomesBinned{x} = Outcomes(find(Confidence == x));
    Mean(x) = nanmean(SubjectOutcomesBinned{x});
    Error(x) = serr(SubjectOutcomesBinned{x});
    ModelOutcomesBinned{x} = so(find(ModelConf_binned == x));
    ModelMean(x) = nanmean(ModelOutcomesBinned{x});
    ModelError(x) = serr(ModelOutcomesBinned{x});
end

% [r2 rmse] = rsquare(Mean,ModelMean);

%% Plot
E2 = errorbar(Mean,Error, 'Color', [.5 .5 .5],'LineWidth',3, 'Marker', '.', 'MarkerSize', 20, 'LineStyle', 'none');
errorbar_tick(E2, 0);
hold on
plot(ModelMean,'LineWidth', 4, 'Color', 'k')
xlabel('Confidence (1-5)', 'FontSize', 24, 'FontName', 'Arial');
set(gca, 'YLim', [.4 1])
set(gca, 'XLim', [0 6])
set(gca, 'XTick', [1 2 3 4 5]);
set(gca, 'XTickLabel', {'1'; '2'; '3'; '4'; '5'})
set(gca, 'YTick', [.5 .75 1]);
set(gca, 'YTickLabel', {'50'; '75'; '100'})
ylabel('% correct', 'FontSize', 24, 'FontName', 'Arial');
set(gca, 'Box', 'off', 'TickDir', 'out', 'FontName', 'Arial', 'TickLength', [0.025 0.05], 'FontSize', 18);
line([0 6], [.5 .5],'Color','k','LineStyle','--', 'LineWidth', 2)



% -------------------------------------------------------------------------
function D = Ddist(N,M,SD)

D = SD * randn(1,N) + M;

function Output = Rail2RailMap(Input)
Output = Input;
Output(Input == 1) = .5;
Output(Input == 2) = .625;
Output(Input == 3) = .75;
Output(Input == 4) = .875;
Output(Input == 5) = 1;

function Output = Map2FivePtScale(Input)
Output(Input < .6) = 1;
Output(Input>.6&Input<.7) = 2;
Output(Input>.7&Input<.8) = 3;
Output(Input>.8&Input<.9) = 4;
Output(Input>.9&Input<1) = 5;