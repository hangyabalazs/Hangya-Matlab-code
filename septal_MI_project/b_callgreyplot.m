function b_callgreyplot
%CALLGREYPLOT   Keypress function for ICA_GUI.
%   Keypress functions for 'Iterative Cluster Analysis' graphical user interface:
%       e - exports the upper figure in a new figure window
%       r - exports the lower figure in a new figure window
%       t - exports the upper accessoric subplot in a new figure window
%       z - exports the lower accessoric subplot in a new figure window
%       y - plots wavelet average magnitudes from 1 to 3 hz against time on the spike train
%       x - plots wavelet average magnitudes from 3 to 6 hz against time on the spike train
%       c - plots wavelet average magnitudes from 6 to 20 hz against time on the spike train
%       v - plots wavelet average magnitudes from 20 to 50 hz against time on the spike train
%       l - draws one line per burst connecting the intraburst action potentials
%
%   See also CALLEXPOTFIG, CALLFFTSWITCH, EXPORTFIG, GREYPLOT and BURSTLINER.

inp = get(gcf,'CurrentCharacter');
switch inp
case 'e'
    b_exportfig('e')
case 'r'
    b_exportfig('r')
case 't'
    b_exportfig('t')
case 'z'
    b_exportfig('z')
case 'y'
    b_greyplot('r')
case 'x'
    b_greyplot('b')
case 'c'
    b_greyplot('g')
case 'v'
    b_greyplot('m')
case 'l'
    b_burstliner
end