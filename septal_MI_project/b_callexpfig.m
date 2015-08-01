function b_callexpfig
%CALLEXPFIG   Keypress function for ICA_GUI.
%   Keypress functions for 'Iterative Cluster Analysis' graphical user interface:
%       e - exports the figure in a new figure window (exports only spike train,
%           return plot or variance plot)
%
%   See also CALLGREYPLOT, CALLFFTSWITCH, CALLEXPORTFIG, EXPFIG and EXPORTFIG.

inp = get(gcf,'CurrentCharacter');
if inp == 'e'
    b_expfig
end