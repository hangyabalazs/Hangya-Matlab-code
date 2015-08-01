function nbpoppsth_call
%NBPOPSTH_CALL   Calls NBPOPPSTH.
%   NBPOPSTH_CALL calls NBPOPPSTH to calculate average PSTHs for identified
%   cholinergic and GABAergic cells.
%
%   See also NBPOPPSTH.

%   Edit log: BH 7/3/12

% I = [825 1171 1176 2051];  % cholinergic
% nbpoppsth(I,false)

% Call NBPOPPSTH for cholinergic cells
I = [1171 1176 2051 1699];   % modified cholinergic
nbpoppsth(I,false)

% Call NBPOPPSTH for GABAergic cells
I = [1455 1456 1511 1512 1515 1583 1655 1680 1681 1716 1764 1859 1885 1967 2165 ...
    1970 2071 ...
    1845 2309];  % GABAergic
nbpoppsth(I,false)