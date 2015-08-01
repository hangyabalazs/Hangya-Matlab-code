function [quociens,dvd] = b_relfreq_for_relfreqrun
%RELFREQ_FOR_RELFREQRUN Version of RELFREQ used by RELFREQRUN.
%   [QUOCIENS,DVD] = RELFREQ_FOR_RELFREQRUN has two output arguments: DVD is diff(vdisc) and
%   QUOCIENS contains the relative frequency.
%
%   See also DISC, DIFF, RELFREQ and RELFREQRUN.

%Input arguments check
error(nargchk(0,0,nargin));

%Input variables
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

%Discrimination variables
global DISC
output = DISC{2};
vdisc = DISC{3};
kuszob = DISC{4};
instfrek = DISC{5};
isi = DISC{6};

%Computing relative frequency
dvd = diff(vdisc);
mn = mean(dvd);
differencia = dvd - mn;
quociens = dvd/mn;

%   List of files calling RELFREQ_FOR_RELFREQRUN:
%   b_relfreqrun