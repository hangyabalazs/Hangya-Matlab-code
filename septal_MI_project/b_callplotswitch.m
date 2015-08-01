function b_callplotswitch
%CALLPLOTSWITCH   Associate keypress function to THETASELECTOR_BETA plots.
%   CALLPLOTSWITCH assigns PLOTSWITCH as keypress function for THETASELECTOR_BETA plots.
%   An easy way to turn PLOTSWITCH on from the Command Window.
%
%   See also PLOTSWITCH and THETASELECTOR_BETA.

set(gcf,'KeyPressFcn','b_plotswitch')