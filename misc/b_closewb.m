function b_closewb
%CLOSEWB   Close progress indicators.
%   CLOSEWB closes all open progress indicators opened by new generation
%   batch processing functions.
%
%   See also CLOSE.

% Get waitbar handles
global WB
wb = WB;

% Close waitbars
for i = 1:length(wb)
    if ishandle(wb(i))
        close(wb(i))
    end
end