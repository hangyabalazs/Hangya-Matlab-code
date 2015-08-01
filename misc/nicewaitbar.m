function WB = nicewaitbar(x,whichbar, varargin)
%NICEWAITBAR   Set waitbar properties.
%   NICEWAITBAR calls waitbar and changes its colors to new color scheme.
%
%   See also WAITBAR.

% Figure properties
WB = waitbar(x,whichbar, varargin{:});   % create waitbar
set(WB,'Color','k')   % black background

% Axes properties
A = allchild(WB);   % axes handle
T = get(A,'title');   % title
set(T,'Color','w')    % title in white

% Patch properties
AA = allchild(A);
aat = get(AA,'Type');
inx = strcmp(aat,'patch');
P = AA(inx);   % patch handle
set(P,'FaceColor',[0.8 0 0])   % patch color