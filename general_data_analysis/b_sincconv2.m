function zint = b_sincconv2(unit,vdisc,newstep)
%SINCCONV2   Convulation with sinc function.
%   Z = SINCCONV2(UNIT,VDISC,NEWSTEP) convulates the input UNIT with sinc function, 
%   when NEWSTEP is the rate between the new and the former sampling frequency, i.e.
%   10 kHz by the resampling of the eeg (it determins the length and the resolution
%   of the output). Discriminated unit (VDISC) is also given, and - unlike SINCCONV -
%   no more thresholding is needed.
%
%   See also SINC and SINCCONV.

% Variable transformation regarding the new sampling frequency
nunit = length(unit);
fs = 10000;       % old sampling frequency in Hz
dt = 1 / fs;
fcut = 100; 
fsnew = 10000 / newstep;
dtnew = newstep / 10000;
fsratio = fsnew / fs;
told = vdisc * dt * fcut;          
tnew = (1:nunit*fsratio) * dtnew * fcut;  

% Convolution
lentold = length(told);
zint = 0;     % final sum
for i = 1:lentold
    zint = zint+sinc(tnew-told(i));
end