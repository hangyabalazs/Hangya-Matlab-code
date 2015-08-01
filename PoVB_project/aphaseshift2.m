function aphaseshift2(chinp,chout)
%APHASESHIFT2   Phaseshift caused by low impedance headstage.
%   APHASESHIFT2(CH1,CH2) calculates the phase difference between 
%   corresponding sine waves in CH1 and CH2. Sine waves are fitted on the 
%   data and phaseshift is calculated using the Hilbert method.
%
%   Syntax:
%   load('Y:\Names\Peter\balazs\phaseshift.mat')
%   aphaseshift2(fazistolasProba1_Ch2,fazistolasProba1_Ch1)
%
%   Reference:
%   Nelson MJ, Pouget P, Nilsen EA, Patten CD, Schall JD (2008) Review of
%   signal distortion through metal microelectrode recording circuits and 
%   filters. J Neurosci Methods 169:141-157. 
%
%   See also APHASESHIFT and SINEFIT.

% Result directory
global DATAPATH
resdir = [DATAPATH 'Andi\Ketxyl\Phaseshift\'];
mm = pwd;

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
    frest = sr / mean(diff(vfchi(k:k+4)));
    ampest = max(linp) - min(linp);
    [xi chisin] = sinefit(linp,sr,[pi/2 frest ampest],[pi/4 frest/2 ampest-1],[pi*3/4 frest*2 ampest+1]);
    ahinp = angle(hilbert(chisin));
    
    lout = fcho(vfchi(k):vfchi(k+4));
    ampest = max(lout) - min(lout);
    [xo chosin] = sinefit(lout,sr,[pi/2 frest ampest],[0 frest/2 ampest-1],[pi*2 frest*2 ampest+1]);
    ahout = angle(hilbert(chosin));
    
    phss = ahinp - ahout;
    frs(next) = xi(2);
    phs(next) = circular_mean(phss,'rad') * 180 / pi;
    next = next + 1;
end

% Plot
figure
plot(frs,phs)
figure
semilogx(frs,phs)

% Save
cd(resdir)
save('phaseshift','frs','phs')
cd(mm)