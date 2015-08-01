function aphaseshift(chinp,chout)
%APHASESHIFT   Phaseshift caused by low impedance headstage.
%   APHASESHIFT(CH1,CH2) calculates the phase difference between 
%   corresponding sine waves in CH1 and CH2. Phaseshift is calculated using
%   the Hilbert method.
%
%   Note, that APHASHIFT2 yields more reliable results!
%
%   Syntax:
%   load('Y:\Names\Peter\balazs\phaseshift.mat')
%   aphaseshift(fazistolasProba1_Ch2,fazistolasProba1_Ch1)
%
%   Reference:
%   Nelson MJ, Pouget P, Nilsen EA, Patten CD, Schall JD (2008) Review of
%   signal distortion through metal microelectrode recording circuits and 
%   filters. J Neurosci Methods 169:141-157. 
%
%   See also APHASESHIFT2.

% Filter
chi = chinp.values;
cho = chout.values;
hlf = round(length(chi)*2/3);
sr = 20000;
nqf = sr / 2;      % Nyquist frequency
flt1 = fir1(4096,40/nqf,'low');      % lowpass filtering below 40 Hz
flt2 = fir1(4096,100/nqf,'low');      % lowpass filtering below 100 Hz
fchi1 = filtfilt(flt1,1,chi(1:hlf));
fcho1 = filtfilt(flt1,1,cho(1:hlf));
fchi2 = filtfilt(flt2,1,chi(hlf+1:end));
fcho2 = filtfilt(flt2,1,cho(hlf+1:end));
fchi = filtfilt(flt2,1,chi);
fcho = filtfilt(flt2,1,cho);

% Standardization
fchi1 = (fchi1 - mean(fchi1)) / std(fchi1);
fchi2 = (fchi2 - mean(fchi2)) / std(fchi2);
fcho1 = (fcho1 - mean(fcho1)) / std(fcho1);
fcho2 = (fcho2 - mean(fcho2)) / std(fcho2);
fchi = (fchi - mean(fchi)) / std(fchi);
fcho = (fcho - mean(fcho)) / std(fcho);

% Discrimination
vfchi1 = disc(-fchi1,0.5);
inx = [];
for k = 1:length(vfchi1)
    if min(fchi1(vfchi1(k)-5000:vfchi1(k)+5000)) < fchi1(vfchi1(k))
        inx = [inx k];
    end
end
vfchi1(inx) = [];
vfchi2 = disc(-fchi2,0.5);
vfchi = [vfchi1 hlf+vfchi2];
figure
plot([fchi1; fchi2])
for k = 1:length(vfchi)
    line([vfchi(k) vfchi(k)],[-5 4],'Color','green')
end

vfcho1 = disc(-fcho1,0.5);
inx = [];
for k = 1:length(vfcho1)
    if min(fcho1(vfcho1(k)-5000:vfcho1(k)+5000)) < fcho1(vfcho1(k))
        inx = [inx k];
    end
end
vfcho1(inx) = [];
vfcho2 = disc(-fcho2,0.5);
vfcho = [vfcho1 hlf+vfcho2];
figure
plot([fcho1; fcho2])
for k = 1:length(vfcho)
    line([vfcho(k) vfcho(k)],[-5 4],'Color','green')
end

% Phaseshift
lenv = length(vfchi);
lle = floor(lenv/5);
frs = zeros(1,lle);
phs = zeros(1,lle);
next = 1;
for k = 1:5:lenv
    linp = fchi(vfchi(k):vfchi(k+4));
    ahee = angle(hilbert(linp));
    phinx = vfcho(k:k+4) - vfchi(k);
    phinx = phinx(phinx>0&phinx<length(ahee));
    if length(phinx) < 3
        keyboard
        error('Technical error 75.')
    end
    phss = ahee(phinx) - pi;
    frs(next) = sr / mean(diff(vfchi(k:k+4)));
    phs(next) = circular_mean(phss,'rad') * 180 / pi;
    next = next + 1;
end
% phs(phs<0) = phs(phs<0) + 360;