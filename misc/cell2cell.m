function C2 = cell2cell(C)
%CELL2CELL   Concatenates cell arrays within a cell.
%   C2 = CELL2CELL(C) accepts a cell array of cell arrays. For each cell of
%   the base array, the inner cells are concatenated.
%
%   See also CELL2MAT.

% Input argument check
error(nargchk(1,1,nargin))
if ~iscell(C)
    C2 = C;
    warning('cell2cell:noncellInput','Non-cell input argument to cell2cell: input returned unchanged.')
    return
end
if ~iscell(C{1})
    C2 = C;
    warning('cell2cell:noncellarrayInput','Input argument to cell2cell is not a cell array of cells: input returned unchanged.')
    return
end

% Conversion
lenC = numel(C);
C2 = cell(size(C));
for k = 1:lenC
    d1 = cellfun(@(s)size(s,1),C{k});   % check dimensions
    d2 = cellfun(@(s)size(s,2),C{k});
    if isempty(d1)    % for empty cells
        C2{k} = C{k};
        continue
    end
    if b_isconstant(d1)
        C2{k} = horzcat(C{k}{1:end});   % concatenation
    elseif b_isconstant(d2)
        C2{k} = vertcat(C{k}{1:end});
    else
        error('Dimension mismatch for cell array concatenation.')
    end
end