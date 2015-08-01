function [out_int,vd] = b_segmentselector(in_int,vdisc)
%SEGMENTSELECTOR    Drops registration segments with less than 100 spikes.
%   [OUT_INT,VD] = SEGMENTSELECTOR(IN_INT_VDISC) requires two input arguments:
%   IN_INT - input segmnets (first and last points after each other) and
%   VDISC - 'vdisc' for the whole registration segment (see DISC for details).
%   It returns OUT_INT - segments with at least 100 spikes and VD - a cell
%   containing 'vdisc' for segments in OUT_INT (indeces relative to the starting
%   point of the whole registration segment).
%
%   See also ICA_NEWAVE, ICA_NEWAVE2, ICA_NEWAVE3, ICA_POWER, ICA_PHASE and 
%   HCN_ICA.

% Input agrument check
error(nargchk(2,2,nargin))

% Main
tooshort = [];
vd = {};
segno = length(in_int) / 2;
for s = 1:segno
    sfsp = in_int(s*2-1);
    slsp = in_int(s*2);
    vvd = vdisc(find(vdisc>sfsp&vdisc<slsp));
    if length(vvd) < 100 %spike number criterium
        tooshort(end+1) = s;
        continue
    end
    vd{end+1} = vvd;
end
in_int(2*tooshort) = 0;
in_int(2*tooshort-1) = 0;
out_int = in_int(find(in_int));
if length(out_int) ~= 2 * length(vd)
    error('356');
end