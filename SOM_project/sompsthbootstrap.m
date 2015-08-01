function [pct1 pct2] = sompsthbootstrap
%SOMPSTHBOOTSTRAP   Bootstrap analysis of PSTH sinilarity.
%   PCT = SOMPSTHBOOTSTRAP calculates bootstrap percentile (PCT) of PV mean
%   similarity to average activated population with respect to a bootstrap
%   sample of all PSTHs aligned to home zone out event.
%
%   [PCT1 PCT2] = SOMPSTHBOOTSTRAP also returns percentile (PCT2) of PV
%   PSTH mean within-group similarity ('homogeneity') as compared to the
%   bootstrap sample.
%
%   See also SOMENTROPY_OPTIMIZE2B.

% Load
% load('C:\Dropbox\KepecsLab\_Sachin\psth normalize\NORMPSTHMATRIX_NONTAGGED_HZOUT.mat')
% load('C:\Dropbox\KepecsLab\_Sachin\psth normalize\NORMPSTHMATRIX_PV_HZOUT.mat')
load('C:\Dropbox\KepecsLab\_Sachin\mPFC\psth similarity\NORMPSTHMATRIX_NONTAGGED_HZIN.mat')
load('C:\Dropbox\KepecsLab\_Sachin\mPFC\psth similarity\NORMPSTHMATRIX_SOM_HZIN.mat')

% PSTH variables
posmod = squeeze(NPSTH_psig(:,:,7:207));
negmod = squeeze(NPSTH_nsig(:,:,7:207));
nomod = squeeze(NPSTH_nonsig(:,:,7:207));
% pv = squeeze(NPSTH_pv(:,:,7:207));
pv = squeeze(NPSTH_som(params.NS_inx,:,7:207));
all = [posmod; negmod; nomod];

% template = mean(posmod);   % template for positive modulation
template = mean(negmod);

% Bootstrap analysis
bsn = 1000;   % number of bootstrap samples
pvn = size(pv,1);   % number of PV cells
bsdist1 = nan(1,bsn);
bsdist2 = nan(1,bsn);
for k = 1:bsn
   rp = randperm(size(all,1));
   bssample = all(rp(1:pvn),:);   % bootstrap sample
   bsdist1(k) = mean(bssample*template');   % bootstrap distribution for similarity w average activated cell
   bsdist2(k) = mean(allsim(bssample));   % bootstrap distribution for average within-group similarity
end

% Get percentile
pvv1 = mean(pv*template');
pct1 = sum(bsdist1>pvv1) / length(bsdist1);
pvv2 = mean(allsim(pv));
pct2 = sum(bsdist2>pvv2) / length(bsdist2);
keyboard

% -------------------------------------------------------------------------
function S = allsim(X)

% Combinations
cb = combnk(1:size(X,1),2);

% Similarities
NumComb = size(cb,1);
S = nan(1,NumComb);
for k = 1:NumComb
    S(k) = dot(X(cb(k,1),:),X(cb(k,2),:));
end