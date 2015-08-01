function newls = copyxls_temp(oldl,newl)
%COPYXLS  Copies information between tables.
%   NEWLS = COPYXLS(OLDL,NEWL) copies information from OLDL to NEWL using
%   the first columns as tags. Columns from matched tags are copied to the
%   new table and returned in NEWLS.
%
%   See alse XLSREAD and XLSWRITE.

% Initialize output
newls = newl;

% Copy based on matching tags
NUMtags = size(oldl,1);
newtags = newl(:,1);  % list of tags in the new table
for iT = 1:NUMtags    % loop through the rows of the old table
    oldtag = oldl{iT,1};   % tag in the old table
    matchedtag = ismember(newtags,oldtag);   % index of the matched tag on the new table
    if ~isempty(oldl{iT,5})
        newls(matchedtag,2) = oldl(iT,5);
    end
    if ~isempty(oldl{iT,6})
        if isequal(oldl{iT,6},2) || isequal(oldl{iT,6},3)
            newls(matchedtag,3) = {0};
        elseif isequal(oldl{iT,6},4)
            newls(matchedtag,3) = {1};
        end
    end
end