function [highfr,instfr] = b_maxf_for_maxfrun
%MAXF_FOR_MAXFRUN   Version of MAXF used by MAXFRUN.
%   [HIGHFR,INSTFR] = MAXF_FOR_MAXFRUN returns the 50 highest value of the instant frequency matrix,
%   excluding values higher than 500 Hz (regarded as artefact) and the instant frequency matrix.
%
%   See also MAXF and MAXFRUN.

% Input arguments check
error(nargchk(0,0,nargin));

% Input variables
global IN
data = IN{1};
eeg = IN{2};
fname = IN{3};
pathname = IN{4};
datinx1 = IN{5};
datinx2 = IN{6};
time = IN{7};
unit = IN{8};
dt = IN{9};
meret = IN{10};
mintafr = IN{11};
xlimit = IN{12};

% Discrimination variables
global DISC
output = DISC{2};
vdisc = DISC{3};
kuszob = DISC{4};
instfrek = DISC{5};
isi = DISC{6};

isi = diff(vdisc);
instfr = (1./isi)*10000;
unreal = find(instfr==max(instfr));
while instfr(unreal) > 500
    instfr(unreal) = 0;
    unreal = find(instfr==max(instfr));
end;
sifr = sort(instfr);
highfr = sifr(end-50:end);