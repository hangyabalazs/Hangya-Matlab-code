function b_callexportfig
%CALLEXPORTFIG   Keypress function for ICA_GUI.
%   Keypress functions for 'Iterative Cluster Analysis' graphical user interface:
%       e - exports the upper figure in a new figure window
%       r - exports the lower figure in a new figure window
%       t - exports the upper accessoric subplot in a new figure window
%       z - exports the lower accessoric subplot in a new figure window
%
%   See also CALLGREYPLOT, CALLFFTSWITCH and EXPORTFIG.

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
end