function b_callhidelines
%CALLHIDELINES   Keypress function for ICA_WAVE contour plots.
%   Keypress functions for ICA_WAVE figures:
%       h - switches the visibility of white lines limiting the wavelet parts for each burst
%           and black frames showing the same wavelet parts excuding the flags before and
%           after bursts.
%
%   See also ICA_WAVE and HIDELINES.

inp = get(gcf,'CurrentCharacter');
if inp == 'h'
    b_hidelines
end