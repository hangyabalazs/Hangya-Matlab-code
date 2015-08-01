function closebut(H)
%CLOSEBUT   Close all figures except specified.
%   CLOSEBUT(H) closes all figures except figure handles H.
%
%   See also CLOSE.

% Get all figure handles
fgs = findobj(allchild(0),'Type','figure');

% Handles to close
clh = setdiff(fgs,H);

% Close
close(clh)