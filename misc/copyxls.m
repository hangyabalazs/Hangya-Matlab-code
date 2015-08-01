function newls = copyxls(oldl,newl)
%COPYXLS  Copies information between tables.
%   NEWLS = COPYXLS(OLDL,NEWL) copies information from OLDL to NEWL using
%   the first columns as tags. Columns from matched tags are copied to the
%   new table and returned in NEWLS.
%
%   See alse XLSREAD and XLSWRITE.

% Initialize output
newls = newl;

% Find out which columns should be copied
cols2add = size(oldl,2) - size(newl,2);
inx1 = size(newl,2) + 1;
inx2 = size(newl,2) + cols2add;

% Copy based on matching tags
NUMtags = size(oldl,1);
newtags = newl(:,1);  % list of tags in the new table
for iT = 1:NUMtags    % loop through the rows of the old table
    oldtag = oldl{iT,1};   % tag in the old table
    matchedtag = ismember(newtags,oldtag);   % index of the matched tag on the new table
    newls(matchedtag,inx1:inx2) = oldl(iT,inx1:inx2);   % copy
end