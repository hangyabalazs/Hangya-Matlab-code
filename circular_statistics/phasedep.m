function [mn_phase mn_rt] = phasedep(spvr,edges)
%PHASEDEP   Phase dependence.
%   [MP MR] = PHASEDEP(S,EDGES) calculates phase dependence of the variable
%   given in S with respect to phase bins determined by EDGES. S should be
%   an N-by-2 matrix with the variable in question in the second column and
%   corresponding phase values in the first column. Mean values conditioned
%   on the phase intervals are returned in MR; sample sizes in the phase
%   bins are returned in MP.
%
%   See also PHASEHIST.

n = length(edges);
mn_phase = zeros(1,n-1);
mn_rt = zeros(1,n-1);
for k = 2:n
    inx = find(spvr(:,1)>edges(k-1)&spvr(:,1)<edges(k));
    mn_phase(k-1) = length(inx);
    mn_rt(k-1) = mean(spvr(inx,2));
end