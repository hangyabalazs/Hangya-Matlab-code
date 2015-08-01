function viptuning
%VIPTUNING   Spectro-temporal tuning curve.
%   VIPTUNING calculates average auditory spectro-temporal tuning curves 
%   for groups of neurons. A scaled Gaussian function is fitted to each
%   smoothed indiviual tuning curve to determine the best frequency. Tuning
%   curves are aligned to their best frequeny and averaged. The average
%   tuning curves are also fitted with Gaussian curves.
%
%   See also VIPPLOTRESPONSE2.

% Load individual tuning curves
load('C:\Balazs\_data\VIP\freq_tuning_data_inh.mat')   % inhibited group
load('C:\Balazs\_data\VIP\freq_tuning_data_act.mat')   % activated group
load('C:\Balazs\_data\VIP\freq_tuning_data_inhact.mat')   % inhibited-activated group
load('C:\Balazs\_data\VIP\freq_tuning_data_VIP.mat')   % VIP neurons
inh = [inh inhact1];

% Average tuning curve for the inhibited group
meantuning(f_tone,inh)

% Average tuning curve for the activated group
meantuning(f_tone,act)

% Average tuning curve for the inhibited-activated group
meantuning(f_tone,inhact1)
meantuning(f_tone,inhact2)

% Average tuning curve for the VIP group
meantuning(f_tone,VIP)

% -------------------------------------------------------------------------
function meantuning(f_tone,inh)

% Fit Gaussian curves
NumTones = length(f_tone);   % number of different tones
NumCells = length(inh);   % number of cells
toneinx = 1:NumTones;   % default index set for tones
[mu_stim, mu_nostim, a_stim a_nostim b_stim b_nostim max_stim, max_nostim] = ...
    deal(nan(1,NumCells));   % preallocate
[fits_stim fits_nostim] = deal(cell(1,NumCells));
[fitted_nostim fitted_stim] = deal({});
for k = 1:NumCells;  % loop throuh all cells 
    
    % Eliminate NaNs
    xc = toneinx;
    nostimrate = inh{k}.rate_tone;   % without laser stimulation
    stimrate = inh{k}.rate_laser;   % with laser stimulation
    naninx = isnan(stimrate) | isnan(nostimrate);
    nostimrate(naninx) = [];   % the fitting programs cannot deal with NaNs
    stimrate(naninx) = [];
    xc(naninx) = [];
    
    % Plot individual tuning curves
    H = figure;
    plot(xc,nostimrate,'ko')
    hold on
    plot(xc,stimrate,'o','MarkerFaceColor','b')

    % Smooth tuning curves
    sstim = smooth(stimrate,'linear',3);  % moving average over 3 neighboring tones
    snostim = smooth(nostimrate,'linear',3);
    
    % Initial mean for the fit
    mloc1 = find(nostimrate==max(nostimrate));
    mloc1 = mloc1(1);   % maximum location for no-stim. tuning curve
    mloc2 = find(stimrate==max(stimrate));
    mloc2 = mloc2(1);   % maximum location for stim. tuning curve
    mloc = mean(mloc1,mloc2);   % use a point in between the two maximum locations
    
    % Fit options
    s = fitoptions('Method','NonlinearLeastSquares',...   % non-linear least squares
        'Lower',[0,0,0,0],...    % lower bounds
        'Upper',[Inf,Inf,xc(end),Inf],...   % upper bounds
        'Startpoint',[max(snostim),0,xc(mloc) 1],...   % initial values
        'Robust','on');    % robust fit
    
    % Model (curve): Gaussian with an additive and a multiplicative constant
    f = getmodel(s);
    
    % Fit
    [fun gof] = fit(xc',snostim,f,s);
    mu_nostim(k) = fun.mu;   % best frequency (no stimulation)
    a_nostim(k) = fun.a;   % scaling parameter
    b_nostim(k) = fun.b;   % constant parameter
    fits_nostim{k} = fun;   % fitted curve
    
    % Plot
    hold on
    plot(xc,snostim,'k:')  % plot tuning curve (no stimulation)
    L = plot(fun,'k');   % overlay the fitted curve
    max_nostim(k) = max(get(L,'YData'));   % maximum of the fitted curve (for normalization)
    fitted_nostim{k} = get(L,'YData');   % fitted curve
    
    % Fit options
    s = fitoptions('Method','NonlinearLeastSquares',...   % non-linear least squares
        'Lower',[0,0,0,0],...   % lower bounds
        'Upper',[Inf,Inf,xc(end),Inf],...   % upper bounds
        'Startpoint',[max(sstim),0,xc(mloc) 1],...   % initial values
        'Robust','on');   % robust fit
    
    % Model (curve): Gaussian with an additive and a multiplicative constant
    f = getmodel(s);
    
    % Fit
    [fun gof] = fit(xc',sstim,f,s);
    mu_stim(k) = fun.mu;   % best frequency (stimulation)
    a_stim(k) = fun.a;   % scaling parameter
    b_stim(k) = fun.b;   % constant parameter
    fits_stim{k} = fun;   % fitted curve
    
    % Plot
    hold on
    plot(xc,sstim,'b:')  % plot tuning curve (stimulation)
    L = plot(fun,'b');   % overlay the fitted curve
    max_stim(k) = max(get(L,'YData'));   % maximum of the fitted curve (for normalization)
    fitted_stim{k} = get(L,'YData');   % fitted curve
    close(H)
end

% Align individual tuning curves to best frequency
[stim_aligned, nostim_aligned] = deal(nan(NumCells,NumTones*2));   % preallocate
xc = toneinx;
for k = 1:NumCells   % loop through all cells
    
    % Interpolate missing frequencies
    nostimrate = inh{k}.rate_tone;   % without laser stimulation
    stimrate = inh{k}.rate_laser;   % with laser stimulation
    nninx = find(~isnan(stimrate)&~isnan(nostimrate));  % interpolate missing frequencies for averaging
    nostimrate = linterp(xc(nninx),nostimrate(nninx),xc);
    stimrate = linterp(xc(nninx),stimrate(nninx),xc);
    nostimrate(nostimrate<0) = 0;  % firing rate should not go below 0
    stimrate(stimrate<0) = 0;   % firing rate should not go below 0
    
    % Align
    lside = sum(xc<mu_nostim(k));   % number of points lower than the best freq.
    rside = NumTones - lside;   % number of points higher than the best freq.
    nostim_aligned(k,(NumTones-lside):(NumTones+rside-1)) = nostimrate;   % align no-stim. tuning curves: move best freq. to the middle
    
    lside = sum(xc<mu_stim(k));   % number of points lower than the best freq.
    rside = NumTones - lside;   % number of points higher than the best freq.
    stim_aligned(k,(NumTones-lside):(NumTones+rside-1)) = stimrate;   % align stim. tuning curves: move best freq. to the middle
end

% Firing rate with and without stimulation
% mn_nostim = nanmean(nostim_aligned,2);
% mn_stim = nanmean(stim_aligned,2);
% mn_nostim = nanmax(nostim_aligned,[],2);  % best option
% mn_stim = nanmax(stim_aligned,[],2);
% mn_nostim = a_nostim;
% mn_stim = a_stim;
% mn_nostim = max_nostim;
% mn_stim = max_stim;
% mn_nostim = nostim_aligned(:,NumTones);
% mn_stim = stim_aligned(:,NumTones);
% figure   % scatter plot for firing rate 'gain'
% plot(mn_nostim,mn_stim,'ko')
% xx = 0:0.1:45;
% yy = xx*0.35;
% hold on
% plot(xx,yy,'r')

% Normalize
% mx = max_nostim';  % normalize to maximum of the fit
mx = max(nostim_aligned,[],2);   % norm. to max. response
% mx = nostim_aligned(:,NumTones);   % norm. to peak response (Scanziani)
nostim_aligned_norm = nostim_aligned ./ repmat(mx,1,size(nostim_aligned,2));   % normalize, no stim.
stim_aligned_norm = stim_aligned ./ repmat(mx,1,size(stim_aligned,2));   % normalize, stim.
ginx = sum(~isnan(stim_aligned_norm))>=5 & sum(~isnan(nostim_aligned_norm))>=5;  % indices where numbers are high enough

% Extrapolate on the flanks
for k = 1:NumCells
    xc2 = 1:NumTones*2;
    nninx = ~isnan(nostim_aligned_norm(k,:));
    nostim_aligned_norm(k,:) = interp1(xc2(nninx),nostim_aligned_norm(k,nninx),xc2,'nearest','extrap');
    nninx = ~isnan(stim_aligned_norm(k,:));
    stim_aligned_norm(k,:) = interp1(xc2(nninx),stim_aligned_norm(k,nninx),xc2,'nearest','extrap');
end

% Average
nostim_mn = nanmean(nostim_aligned_norm);   % average tuning curve, no stim.
stim_mn = nanmean(stim_aligned_norm);   % average tuning curve, stim.
nostim_err = nanstd(nostim_aligned_norm) ./ sqrt(sum(~isnan(nostim_aligned_norm)));   % SE, no stim.
stim_err = nanstd(stim_aligned_norm) ./ sqrt(sum(~isnan(stim_aligned_norm)));   % SE, stim.

% Eliminate points with less than 5 numbers in the average
xc = 1:length(nostim_mn);   % for x axis
naninx = ~ginx;  % indices to eliminate
nostim_mn(naninx) = [];
stim_mn(naninx) = [];
nostim_err(naninx) = [];
stim_err(naninx) = [];
xc(naninx) = [];

% Initial mean for the fit
mloc1 = find(nostim_mn==max(nostim_mn));
mloc1 = mloc1(1);   % maximum location for average no-stim. tuning curve
mloc2 = find(stim_mn==max(stim_mn));
mloc2 = mloc2(1);   % maximum location for average stim. tuning curve
mloc = mean(mloc1,mloc2);   % use a point in between the two maximum locations

% Fit options
s = fitoptions('Method','NonlinearLeastSquares',...   % non-linear least squares
    'Lower',[0,0,0,0],...   % lower bounds
    'Upper',[Inf,Inf,xc(end),Inf],...   % upper bounds
    'Startpoint',[max(nostim_mn),0,xc(mloc) 1],...   % initial values
    'Robust','on');   % robust fit

% Model (curve): Gaussian with an additive and a multiplicative constant
f = getmodel(s);

% Fit
[fun_nostim gof] = fit(xc',nostim_mn',f,s);

% Fit options
s = fitoptions('Method','NonlinearLeastSquares',...   % non-linear least squares
    'Lower',[0,0,0,0],...  % lower bounds
    'Upper',[Inf,Inf,xc(end),Inf],...   % upper bounds
    'Startpoint',[max(stim_mn),0,xc(mloc) 1],...   % initial values
    'Robust','on');  % robust fit

% Model (curve): Gaussian with an additive and a multiplicative constant
f = getmodel(s);

% Fit
[fun_stim gof] = fit(xc',stim_mn',f,s);

% Plot
figure
E1 = errorbar(xc,nostim_mn,nostim_err,'k.','MarkerSize',16,'LineWidth',2);
hold on
E2 = errorbar(xc,stim_mn,stim_err,'b.','MarkerSize',16,'LineWidth',2);
L1 = plot(fun_nostim,'k');
L2 = plot(fun_stim,'b');
errorbar_tick(E1,0)
errorbar_tick(E2,0)
set([L1 L2],'LineWidth',2)
xlim([xc(1)-0.5 xc(end)+0.5])
box off

keyboard

% Overlay scaled/shifted version of the control tuning curve
% Fit options
s = fitoptions('Method','NonlinearLeastSquares',...   % non-linear least squares
    'Lower',-Inf,...  % lower bound
    'Upper',Inf,...   % upper bound
    'Startpoint',1,...   % initial value
    'Robust','on');  % robust fit

% Model (curve): scaling by a multiplicative constant
f = fittype(['a*(' num2str(fun_nostim.a) '*(exp(-(x-' ...
    num2str(fun_nostim.mu) ')^2/(2*' num2str(fun_nostim.sigma) '^2)))+' ...
    num2str(fun_nostim.b) ')'],'options',s);

% Fit
[fun_stim_scale gof_scale] = fit(xc',stim_mn',f,s);
L3 = plot(fun_stim_scale,'r');
set(L3,'LineWidth',2)

% Fit options
s = fitoptions('Method','NonlinearLeastSquares',...   % non-linear least squares
    'Lower',-Inf,...  % lower bound
    'Upper',Inf,...   % upper bound
    'Startpoint',0,...   % initial value
    'Robust','on');  % robust fit

% Model (curve): shifting by an additive constant
f = fittype(['a+(' num2str(fun_nostim.a) '*(exp(-(x-' ...
    num2str(fun_nostim.mu) ')^2/(2*' num2str(fun_nostim.sigma) '^2)))+' ...
    num2str(fun_nostim.b) ')'],'options',s);

% Fit
[fun_stim_shift gof_shift] = fit(xc',stim_mn',f,s);
L4 = plot(fun_stim_shift,'c');
set(L4,'LineWidth',2)

% Fit options
s = fitoptions('Method','NonlinearLeastSquares',...   % non-linear least squares
    'Lower',[-Inf -Inf],...  % lower bound
    'Upper',[Inf Inf],...   % upper bound
    'Startpoint',[1 0],...   % initial value
    'Robust','on');  % robust fit

% Model (curve): scaling by a multiplicative constant
f = fittype(['a*(' num2str(fun_nostim.a) '*(exp(-(x-' ...
    num2str(fun_nostim.mu) ')^2/(2*' num2str(fun_nostim.sigma) '^2)))+' ...
    num2str(fun_nostim.b) ')+b'],'options',s);

% Fit
[fun_stim_combined gof_combined] = fit(xc',stim_mn',f,s);
L5 = plot(fun_stim_combined,'m');
set(L5,'LineWidth',2)

keyboard

% Model selection for the scale model
n = length(stim_mn);   % number of data points
k = 1;   % number of model parameters
R = gof_scale.rmse .^ 2;   % mean square error
AIC_scale = n * log(R/n) + 2 * k;   % Akaike information criterion
AICc_scale = AIC_scale + (2 * k * (k + 1)) / (n - k - 1);   % AIC corrected for finite sample sizes
BIC_scale = n * log(R/n) + k * log(n);   % Bayesian information criterion

% Model selection for the shift model
n = length(stim_mn);   % number of data points
k = 1;   % number of model parameters
R = gof_shift.rmse .^ 2;   % mean square error
AIC_shift = n * log(R/n) + 2 * k;   % Akaike information criterion
AICc_shift = AIC_shift + (2 * k * (k + 1)) / (n - k - 1);   % AIC corrected for finite sample sizes
BIC_shift = n * log(R/n) + k * log(n);   % Bayesian information criterion

% Model selection for the combined model
n = length(stim_mn);   % number of data points
k = 2;   % number of model parameters
R = gof_combined.rmse .^ 2;   % mean square error
AIC_combined = n * log(R/n) + 2 * k;   % Akaike information criterion
AICc_combined = AIC_combined + (2 * k * (k + 1)) / (n - k - 1);   % AIC corrected for finite sample sizes
BIC_combined = n * log(R/n) + k * log(n);   % Bayesian information criterion

keyboard

% -------------------------------------------------------------------------
function f = getmodel(s)

f = fittype('a*(exp(-(x-mu)^2/(2*sigma^2)))+b','options',s);