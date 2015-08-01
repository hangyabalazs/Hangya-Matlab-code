function maximize_figure(fig)
%MAXIMIZE_FIGURE   Maximize figure window.
%   MAXIMIZE_FIGURE(H) maximizes the figure with handle H. MAXIMIZE_FIGURE
%   with no input arguments maximizes the current figure.
%
%   See also FIGURE and GCF.

%   Copyright to Peter Bartho.

% Use current figure if no input arguments
if nargin == 0
    fig = gcf;
end

% Maximize figure window
units = get(fig,'units');
set(fig,'units','normalized','outerposition',[0 0 1 1]);
set(fig,'units',units);