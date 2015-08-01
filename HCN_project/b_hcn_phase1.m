function b_hcn_phase1
%HCN_PHASE1   Cumulative histogram.
%   Calculates cumulative histogram of HCN analysis phase results without any
%   normalization. Result is plotted (without saving).
%
%   Note: normalization would be necessaray!
%
%   See also HCN_ANALYSIS.

% Input arguments check
error(nargchk(0,0,nargin))

% Input directory
cd d:\_analysis\matlab_data\HCN\Analysis4

% Cumulative histogram
d = dir;
edges = [-pi : pi/30 : pi];
edges = edges(2:end);
AH = zeros(1,length(edges));
for i = 3:length(d)
    if ~d(i).isdir
        if ~isempty(findstr(d(i).name,'PHASEDIFF'))
            load(d(i).name)
            AH = AH + AngleHistogram(1:end-1);
        end
    end
end
bar(edges,AH);