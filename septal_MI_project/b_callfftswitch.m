function b_callfftswitch
%CALLFFTSWITCH   Keypress function for ICA_GUI.
%   Keypress functions for 'Iterative Cluster Analysis' graphical user interface:
%       e - exports the upper figure in a new figure window
%       r - exports the lower figure in a new figure window
%       t - exports the upper accessoric subplot in a new figure window
%       z - exports the lower accessoric subplot in a new figure window
%       q - switches up one wavelet average magnitude fft plot
%       w - switches down one wavelet average magnitude fft plot
%       z - switches 'zoom xon' on
%
%   See also CALLEXPORTFIG, CALLGREYPLOT, EXPORTFIG and FFTSWITCH.

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
case 'q'
    b_fftswitch('u')
case 'w'
    b_fftswitch('d')
case 'z'
    zoom xon
end