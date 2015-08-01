function b_callexportfig0
%CALLEXPORTFIG   Keypress function for ICA_GUI.
%   Keypress functions for 'Iterative Cluster Analysis' graphical user interface:
%       e - exports the upper figure in a new figure window
%       r - exports the lower figure in a new figure window
%   See also CALLGREYPLOT, CALLFFTSWITCH and EXPORTFIG.

inp = get(gcf,'CurrentCharacter');
switch inp
case 'e'
    b_exportfig('e')
case 'r'
    b_exportfig('r')
end