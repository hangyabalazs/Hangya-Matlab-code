function modelcorr(lim1,lim2)
%COND_ACCURACY_FR2   Accuracy conditioned on firing rate.
%   COND_ACCURACY_FR2 calculates the percentage of hits relative to the sum
%   of hits and false alarms for specific firing rate intervals of a given
%   cell. Firing rates are binned; hits and false alarms are partitioned
%   according to the firing rate bins. It does not handle multiple
%   sessions; different stimulus intensities are pooled together. To a more
%   sophisticated, however slower solution, see COND_ACCURACY_FR.
%
%   CONDPERF = COND_ACCURACY_FR2(LIM1,LIM2) uses LIM1 and LIM2 to determine
%   the time window before stimuli used for firing rate calculation (in
%   seconds). Default is -2:0. Conditional performance vector (CONDPERF) is
%   returned.
%
%   See also COND_ACCURACY_FR.

%   Edit log: BH 7/7/11

% Input variables
if nargin == 0
    lim1 = 0;
    lim2 = 0.2;
end
cellid = 'n023_111218a_1.2';

% Load trial events
EventName = 'LeftPortIn';
ST = loadcb(cellid,'EVENTSPIKES');   % load prealigned spikes for stimulation events
TE = loadcb(cellid,'TrialEvents');
event_pos = findcellstr(ST.events(:,1),EventName);
spikes_stimon = ST.event_stimes{event_pos};

% Load model parameters
load('C:\Balazs\_analysis\NB\bimodalITI\fit\model2_ML_nb023_allbimodal_smns_intensities_20_30_40.mat')

% Poststimulus frequency
NUMtrials = length(spikes_stimon);
prestimfreq = nan(1,NUMtrials);
for k = 1:NUMtrials-1
    lspikes = spikes_stimon{k};
    lspikes2 = lspikes(lspikes>lim1&lspikes<lim2);   % time window: one sec before stimulus onset
    prestimfreq(k) = length(lspikes2) / (lim2 - lim1);
end

% Expectations based on the hazard rate model
itis = TE.ITIDistribution;      % foreperiods
itis = itis(1:end-1);
fitexp = lmodel(param,itis);    % temporal expectation based on the fitted model

% Correlation
hitinx = find(~isnan(TE.Hit));  % hits
hfitexp = fitexp(hitinx);
hfr = prestimfreq(hitinx);

H = figure;
plot(hfitexp,hfr,'.')
[b,bint,r,rint,stats] = regress(hfr',[ones(length(hfitexp),1),hfitexp']);
R = sqrt(stats(1));         % correlation coefficient (R-value of the regression)
F = stats(2);           % F-test for H0: all coeff.-s are zero
p = stats(3);           % F-test significance

keyboard

% -------------------------------------------------------------------------
function y = lmodel(param,x)
%LMODEL   Bimodal subjective hazard rate model.
%   LMODEL is modified from BIMODAL_HAZARD_MODEL2; it only returns the
%   values of the bimodal anticipation function.
%       y(t) = Ab(t)
%   where the bimodal anticipation function denoted by 'Ab'. Model
%   parameters should be passed to the function in the PARAM structure. The
%   bimodal anticipation function (Ab) is a mixture of two Gaussians
%   (means, 0.3 and 2s; SD, 0.15s) and a uniform random variable
%   constrained between 0.1 and 3 seconds (mixing proportions: 0.35, 0.35,
%   0.3). Subjective hazard rate is calculated from the hazard rate by
%   convolving it with a Gaussian pdf with variability linearly increasing
%   in time (Phi, slope parameter of this increase).
%
%   See also BIMODAL_HAZARD_MODEL2.

% Handle both vector and structure input
if ~isstruct(param)
    st.we = param(1);
    st.wr = param(2);
    st.wb = param(3);
    st.phi = param(4);
    st.tau = param(5);
    param = st;
end

% Call the model
y = subjective_hazard_rate(param.phi,x-param.tau);

% -------------------------------------------------------------------------
function v = subjective_hazard_rate(Phi,iti)

% Parameters of the foreperiod distribution
ITIMin = 0.1;   % parameters of the uniform distribution
ITIMax = 3;
mng1 = 0.3;   % parameters for the Gaussians
mng2 = 2;
sdg = 0.15;
pmx1 = 0.35;   % mixing probabilities
pmx2 = 0.35;
pmx3 = 1 - pmx1 - pmx2;

% Time
dt = 0.01;     % resolution
t = ITIMin:dt:ITIMax;   % time

% Failure density function
ft = pmx3 / (ITIMax - ITIMin) + ...
    pmx1 * normpdf(t,mng1,sdg) + ...
    pmx2 * normpdf(t,mng2,sdg);
% figure
% plot(t,ft)
ft = ft / sum(ft);  % for estimating cdf, we need a histogram estimation

% Hazard rate
Ft = cumsum(ft);
Rt = 1 - Ft;
ht = (Rt(1:end-1) - Rt(2:end)) ./ (dt * Rt(1:end-1));
% figure
% plot(t(1:end-1),ht)

% Subjective hazard rate:
% Convolution of the failure density function with a time-dependent Gaussian
intfg = zeros(size(t));
tau = t;
cntr = 1;   % counter
t2 = [fliplr(t(1)-dt:-dt:dt) t t(end)+dt:dt:3*t(end)];   % we need a longer time range, to get a better estimate of Fwt
for tm = t2
    fg = ft .* exp(-(tau-tm).^2/(2*(Phi*tm).^2));
    intfg(cntr) = sum(fg*dt);
    cntr = cntr + 1;
end
fwt = 1 ./ (Phi * t2 * sqrt(2*pi)) .* intfg;
fwt = fwt / sum(fwt);

% Subjective hazard rate
Fwt = cumsum(fwt);
Rwt = 1 - Fwt;
hwt = (Rwt(1:end-1) - Rwt(2:end)) ./ (dt * Rwt(1:end-1));
% figure
% plot(t2(1:end-1),hwt)

% Values at requested times
v = linterp(t2(1:end-1),hwt,iti);