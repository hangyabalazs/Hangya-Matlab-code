function b_relfreq
%RELFREQ Calculates relative interspike interval length.
%   The function draws two plots: the one titled 'differentia' is the difference,
%   the other one (titled 'quociens') is the ratio of the interval length and
%   the mean of them.
%
%   See also RELFREQRUN.

%Input arguments check
error(nargchk(0,0,nargin));

%Import
b_in3
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

%Discrimination
b_disc;
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

%Plotting
close all;
figure;
subplot(2,1,1);
plot(differencia);
subplot(2,1,2);
plot(unit);
title('differencia');
figure;
subplot(2,1,1);
plot(quociens);
subplot(2,1,2);
plot(unit);
title('quociens');

%   List of files calling RELFREQ: