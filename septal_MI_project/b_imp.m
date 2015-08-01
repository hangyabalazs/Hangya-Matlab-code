function b_imp(fname,where,data,datinx1,datinx2)
%IMP    Creates import variables within functions.
%   IMP(FNAME,WHERE,DATA,DATINX1,DATINX2) creates the variables 'unit', 'eeg', 'mintafr'
%   (sampling rate), 'dt' (inverse of sampling rate), 'time' (domain of data), 'xlimit'
%   (time limits), 'meret' (size of data) in the workspace of the caller function. It
%   also assignes the global variable IN. The input arguments are FNAME (file name),
%   WHERE (path name), DATA (loaded data), DATINX1 and DATINX2 (limits of the input
%   interval)
%
%   See also IN3 and VAR2WS.

% Input arguments check
error(nargchk(5,5,nargin));

% Creating the import variables
unit = data(datinx1:datinx2,2);
unit = unit';
eeg = data(datinx1:datinx2,1);
eeg = eeg';
dt = 0.0001;
time = [0:length(unit)-1] * dt;
xlimit = [min(time),max(time)];
meret = size(data,1);
mintafr = 1 / dt;

% Create global IN
global IN
IN = cell(1,12);
IN{1} = data;
IN{2} = eeg;
IN{3} = fname;
IN{4} = where;   %path name
IN{5} = datinx1;
IN{6} = datinx2;
IN{7} = time;
IN{8} = unit;
IN{9} = dt;
IN{10} = meret;
IN{11} = mintafr;
IN{12} = xlimit;

% Assigning the variables in the caller workspace
ws = 'caller';
assignin(ws,'eeg',IN{2})
assignin(ws,'time',IN{7})
assignin(ws,'unit',IN{8})
assignin(ws,'dt',IN{9})
assignin(ws,'meret',IN{10})
assignin(ws,'mintafr',IN{11})
assignin(ws,'xlimit',IN{12})