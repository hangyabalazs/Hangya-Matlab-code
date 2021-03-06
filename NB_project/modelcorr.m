function [R p H] = modelcorr(cellid,varargin)
%MODELCORR   Correlation of firing rate and anticipation.
%   [R P H] = MODELCORR(CELLID) calculates linear correlation (R;
%   significance, P) between post-stimulus firing rate (window, 0-0.5 s)
%   and temporal expectancy for a given cell (CELLID). Temporal expectancy
%   is estimated based on a hazard rate model fitted to the behavior;
%   fitting is performed outside this program and only the model parameters
%   are loaded here (see BIMODAL_HAZARD_MODEL2). Scatter plot with line fit
%   is drawn in figure H ('display' has to be set to true, see below).
%   Default behavior can be overwritten using the following parameter-value
%   pairs as optional input arguments (with default values):
%       'window',[0 0.5] - timing relative to the event for spike rate
%           window
%       'event', 'LeftPortIn' - reference event for spike rate window
%       'display', false - control plotting behavior
%
%   See also RTCORR.

% Default arguments
prs = inputParser;
addRequired(prs,'cellid',@iscellid)
addParamValue(prs,'window',[0 0.5],@(s)isnumeric(s)&isequal(length(s),2))  % time window relative to the event, in seconds
addParamValue(prs,'event','LeftPortIn',@ischar)   % default reference event: 'LeftPortIn'
addParamValue(prs,'display',false,@(s)islogical(s)|ismember(s,[0 1]))   % control displaying rasters and PSTHs
parse(prs,cellid,varargin{:})
g = prs.Results;

% Load trial events
ST = loadcb(cellid,'EVENTSPIKES');   % load prealigned spikes for stimulation events
TE = loadcb(cellid,'TrialEvents');
event_pos = findcellstr(ST.events(:,1),g.event);
spikes_stimon = ST.event_stimes{event_pos};

% Load model parameters
animalname = cellid2tags(cellid);
global DATAPATH
dr = [DATAPATH '\NB\bimodalITI\fit\'];
ddr = dir(dr);
ddrf = {ddr.name};
fn = ['model2_ML_' animalname '_allbimodal_smns'];
fninx = find(cellfun(@(s)~isempty(s),strfind(ddrf,fn)));
fnm = [ddrf{fninx(1)}(1:end-4) '.mat'];
load([dr fnm])

% Poststimulus frequency
lim1 = g.window(1);
lim2 = g.window(2);
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

% Filter trials
hitinx = find(~isnan(TE.Hit));  % hits
hfitexp = fitexp(hitinx);
hfr = prestimfreq(hitinx);

% Plot
x = hfitexp;   % temporal expectation
y = hfr;
if g.display
    H = figure;
    plot(x,y,'ko','MarkerFaceColor','k')
    % [gr icp err] = linefit(x,y);
    coef = polyfit(x,y,1);
    gr = coef(1);
    icp = coef(2);
    xx = (min(x)-0.1):0.01:(max(x)+0.1);
    yy = xx .* gr + icp;
    hold on
    plot(xx,yy,'LineWidth',2,'Color','black')
    xlabel('Reaction time (ms)')
    ylabel('Firing rate (Hz)')
else
    H = NaN;
end

% Correlation
[b,bint,r,rint,stats] = regress(hfr',[ones(length(hfitexp),1),hfitexp']); %#ok<*ASGLU>
pR = corrcoef(hfr,hrt);
R = pR(2);
% R2 = sqrt(stats(1));         % correlation coefficient (R-value of the regression)
F = stats(2);           %#ok<NASGU> % F-test for H0: all coeff.-s are zero
p = stats(3);           % F-test significance
disp([R p])

% keyboard

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