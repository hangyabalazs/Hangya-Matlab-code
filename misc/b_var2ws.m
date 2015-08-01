function varargout = b_var2ws(X,ws)
%VAR2WS Loads IN and DISC global variables to workspace.
%   VAR2WS loads global variables DATADIR, DATAPATH, IN and DISC to workspace, 
%   if they exist. It also obtains the data from cell arrays IN and DISC.
%
%   VAR2WS(X) loads only DATADIR, DATAPATH and X, where X should be 'in', 'disc'
%   or 'all'. In case X is 'all', both IN and DISC are loaded.
%
%   VAR2WS(X,WS) loads the variables to the workspace WS, where WS should be
%   either 'base' or 'caller'.
%
%   See also IN, DISC and STARTUP.

% Input arguments check
error(nargchk(0,2,nargin));

if nargin > 0
    switch X
    case 'in'
        in = 1;
        disc = 0;
    case 'disc'
        in = 0;
        disc = 1;
    case 'all'
        in = 1;
        disc = 1;
    otherwise
        error('Input must be ''in'' or ''disc''.');
    end
else
    in = 1;
    disc = 1;
end

if nargin < 2
    ws = 'base';
end

% Loading IN
if in
    global IN
    if isempty(IN)
        clear global IN
    else
        assignin(ws,'data',IN{1})
        assignin(ws,'eeg',IN{2})
        assignin(ws,'fname',IN{3})
        assignin(ws,'pathname',IN{4})
        assignin(ws,'datinx1',IN{5})
        assignin(ws,'datinx2',IN{6})
        assignin(ws,'time',IN{7})
        assignin(ws,'unit',IN{8})
        assignin(ws,'dt',IN{9})
        assignin(ws,'meret',IN{10})
        assignin(ws,'mintafr',IN{11})
        assignin(ws,'xlimit',IN{12})
    end
end

% Loading DISC
if disc
    global DISC
    if isempty(DISC)
        clear global DISC
    else
        assignin(ws,'id',DISC{1})
        assignin(ws,'output',DISC{2})
        assignin(ws,'vdisc',DISC{3})
        assignin(ws,'kuszob',DISC{4})
        assignin(ws,'instfrek',DISC{5})
        assignin(ws,'isi',DISC{6})
    end
end

% Loading DATADIR
global DATADIR
if isempty(DATADIR)
    clear global DATADIR
else
    assignin(ws,'DATADIR',DATADIR)
end

% Loading DATAPATH
global DATAPATH
if isempty(DATAPATH)
    clear global DATAPATH
else
    assignin(ws,'DATAPATH',DATAPATH)
end