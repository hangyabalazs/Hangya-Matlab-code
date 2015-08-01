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
load('C:\Balazs\_data\VIP\freq_tuning_data_VIP.mat')   % VIP neurons

% Average tuning curve for the inhibited group
meantuning(f_tone,inh)

% -------------------------------------------------------------------------
function meantuning(f_tone,inh)

% Fit Gaussian curves
NumTones = length(f_tone);   % number of different tones
NumCells = length(inh);   % number of cells
toneinx = 1:NumTones;   % default index set for tones
[mu_stim, mu_nostim, max_stim, max_nostim] = deal(nan(1,NumCells));   % preallocate
[fitted_nostim fitted_stim] = deal({});
for k = 1:NumCells;  % loop throuh all cells 
    
    % Eliminate NaNs
    xc = toneinx;
    nostimrate = inh{k}.rate_tone;   % without laser stimulation
    stimrate = inh{k}.rate_laser(:,k);   % with laser stimulation
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
    
    % Plot
    hold on
    plot(xc,sstim,'b:')  % plot tuning curve (stimulation)
    L = plot(fun,'b');   % overlay the fitted curve
    max_stim(k) = max(get(L,'YData'));   % maximum of the fitted curve (for normalization)
    fitted_stim{k} = get(L,'YData');   % fitted curve
    close(H)
end

% Align individual tuning curves to best frequency
NumTones = length(fitted_nostim{1});
[stim_aligned, nostim_aligned] = deal(nan(NumCells,NumTones*2));   % preallocate
for k = 1:NumCells   % loop through all cells
    
    % Interpolate missing frequencies
    nostimrate = fitted_nostim{k};   % without laser stimulation
    stimrate = fitted_stim{k};   % with laser stimulation
    nninx = find(~isnan(stimrate)&~isnan(nostimrate));  % interpolate missing frequencies for averaging
    nostimrate = linterp(xc(nninx),nostimrate(nninx),xc);
    stimrate = linterp(xc(nninx),stimrate(nninx),xc);
    
    % Align
    mloc = find(nostimrate==max(nostimrate));  % find peak
    lside = sum(xc<mloc);   % number of points lower than the best freq.
    rside = NumTones - lside;   % number of points higher than the best freq.
    nostim_aligned(k,(NumTones-lside):(NumTones+rside-1)) = nostimrate;   % align no-stim. tuning curves: move best freq. to the middle
    
    mloc = find(stimrate==max(stimrate));  % find peak
    lside = sum(xc<mloc);   % number of points lower than the best freq.
    rside = NumTones - lside;   % number of points higher than the best freq.
    stim_aligned(k,(NumTones-lside):(NumTones+rside-1)) = stimrate;   % align stim. tuning curves: move best freq. to the middle
end

% Normalize
mx = max_nostim';  % normalize to maximum of the fit
% mx = max(nostim_aligned,[],2);   % norm. to max. response
% mx = nostim_aligned(:,NumTones);   % norm. to peak response (Scanziani)
nostim_aligned_norm = nostim_aligned ./ repmat(mx,1,size(nostim_aligned,2));   % normalize, no stim.
stim_aligned_norm = stim_aligned ./ repmat(mx,1,size(stim_aligned,2));   % normalize, stim.
nostim_mn = nanmean(nostim_aligned_norm);   % average tuning curve, no stim.
stim_mn = nanmean(stim_aligned_norm);   % average tuning curve, stim.
nostim_err = nanstd(nostim_aligned_norm) ./ sqrt(sum(~isnan(nostim_aligned_norm)));   % SE, no stim.
stim_err = nanstd(stim_aligned_norm) ./ sqrt(sum(~isnan(stim_aligned_norm)));   % SE, stim.

% Eliminate points with less than 5 numbers in the average
xc = 1:length(nostim_mn);   % for x axis
ginx = sum(~isnan(stim_aligned_norm))>=5 & sum(~isnan(nostim_aligned_norm))>=5;  % indices where numbers are high enough
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
errorbar(xc,nostim_mn,nostim_err,'ko')
hold on
errorbar(xc,stim_mn,stim_err,'bo')
L = plot(fun_nostim,'k');
plot(fun_stim,'b')

keyboard

xd=get(L,'XData');
yd=get(L,'YData');
plot(xd,yd*0.4,'r')
plot(xd,yd*0.35,'r')

% -------------------------------------------------------------------------
function f = getmodel(s)

f = fittype('a*(exp(-(x-mu)^2/(2*sigma^2)))+b','options',s);