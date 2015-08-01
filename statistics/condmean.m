function [cond_mean cond_median cond_samplesize all_dep all_indep] = ...
    condmean(spvr,edges)
%CONDMEAN   Conditional mean.
%   CM = CONDMEAN(S,EDGES) calculates the dependence of the variable S(:,2)
%   with respect to S(:,1) using bins determined by EDGES. Thus, S should
%   be an N-by-2 matrix with the dependent variable in the second column
%   and the corresponding independent values in the first column. Mean
%   values conditioned on the independent variable are returned in CM.
%
%   [CM1 CM2 CS] = CONDMEAN(S,EDGES) returns conditional mean (CM1), median
%   (CM2) and sample sizes of the bins (CS). NANMEAN and NANMEDIAN are used
%   to calculate the statistics.
%
%   [CM1 CM2 CS AD AI] = CONDMEAN(S,EDGES) also returns dependent (AD) and
%   independent (AI) values sorted into the bins.
%
%   CONDMEAN also supports overlapping bins. In this case, EDGES should be
%   a 2-by-N matrix with corresponding bin limits in the columns.
%
%   See also PHASEDEP3.

% Decide whether overlapping bins are used
[ne me] = size(edges);
if ne > me
    edges = edges';
    [ne me] = size(edges);
end

% Calculate conditional mean values
switch ne
    case 1      % non-overlapping bins
        n = length(edges);
        cond_samplesize = nan(1,n-1);
        cond_mean = nan(1,n-1);
        cond_median = nan(1,n-1);
        all_dep = cell(1,n-1);
        all_indep = cell(1,n-1);
        for k = 2:n
            inx = find(spvr(:,1)>edges(k-1)&spvr(:,1)<=edges(k));
            cond_samplesize(k-1) = length(inx);     % sample size of the bin
            cond_mean(k-1) = nanmean(spvr(inx,2));     % mean in the bin
            cond_median(k-1) = nanmedian(spvr(inx,2)); % median in the bin
            all_dep{k-1} = spvr(inx,2);             % all dependent values in the bin
            all_indep{k-1} = spvr(inx,1);           % corresponding independent values
        end
    case 2      % overlapping bins
        n = me;
        cond_samplesize = nan(1,n);
        cond_mean = nan(1,n);
        cond_median = nan(1,n);
        all_dep = cell(1,n);
        all_indep = cell(1,n);
        for k = 1:n
            inx = find(spvr(:,1)>edges(1,k)&spvr(:,1)<=edges(2,k));
            cond_samplesize(k) = length(inx);     % sample size of the bin
            cond_mean(k) = nanmean(spvr(inx,2));     % mean in the bin
            cond_median(k) = nanmedian(spvr(inx,2)); % median in the bin
            all_dep{k} = spvr(inx,2);             % all dependent values in the bin
            all_indep{k} = spvr(inx,1);           % corresponding independent values
        end
    otherwise
        error('Second input argument must be 1-by-N or 2-by-N.')
end