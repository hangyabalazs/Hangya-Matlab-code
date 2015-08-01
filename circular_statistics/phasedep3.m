function [mn_phase mn_rt] = phasedep3(spvr,edges)
%PHASEDEP3   Phase dependence.
%   [MP MR] = PHASEDEP3(S,EDGES) calculates phase dependence of the variable
%   given in S with respect to phase bins determined by EDGES. S should be
%   an N-by-2 matrix with the variable in question in the second column and
%   corresponding phase values in the first column. Mean values conditioned
%   on the phase intervals are returned in MR; sample sizes in the phase
%   bins are returned in MP.
%
%   PHASEDEP3 also supports overlapping bins. In this case, EDGES should be
%   a 2-by-N matrix with corresponding bin limits in the columns.
%
%   See also PHASEHIST and PHASEDEP.

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
        mn_phase = zeros(1,n-1);
        mn_rt = zeros(1,n-1);
        for k = 2:n
            inx = find(spvr(:,1)>edges(k-1)&spvr(:,1)<edges(k));
            mn_phase(k-1) = length(inx);
            mn_rt(k-1) = mean(spvr(inx,2));
        end
    case 2      % overlapping bins
        n = me;
        mn_phase = zeros(1,n);
        mn_rt = zeros(1,n);
        for k = 1:n
            if edges(2,k) <= pi
                inx = find(spvr(:,1)>edges(1,k)&spvr(:,1)<edges(2,k));
            else    % circular case for edges
                inx = find((spvr(:,1)>edges(1,k)&spvr(:,1)<edges(2,k))|... 
                (spvr(:,1)>(edges(1,k)-2*pi)&spvr(:,1)<(edges(2,k)-2*pi)));
            end
            mn_phase(k) = length(inx);
            mn_rt(k) = mean(spvr(inx,2));
        end
    otherwise
        error('Second input argument must be 1-by-N or 2-by-N.')
end