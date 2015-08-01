function b_rescalefig(scale);
%RESCALEFIG   Changes ticklabeles from 'scale' to 'frequency'.
%   RESCALEFIG(SCALE) requires an input parameter SCALE containing the correspondig
%   frequency value for each scale index.
%
%   See also RETIMEFIG and WAVELET.

a0 = get(gca,'yticklabel');
a1 = scale(str2num(a0));
set(gca,'yticklabel',num2str(a1'))