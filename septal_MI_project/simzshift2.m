function zmaxloc = simzshift2(sr,len,np,tp,delay,thetafreq,uf,unr)
%SIMZSHIFT   Simulations for testing the validity of Z-shift method.
%   SIMZSHIFT calculates Z-shift on simualted artifical unit and LFP data.
%   A random Poisson point process is generated to simulate single cell
%   firing. Gaussian noise is added to a sine wave to mimic a periodic LFP
%   oscillation. Information transmission is established by adding the unit
%   convolved by a Blackman-Harris window to the LFP after applying a given
%   delay and scale parameter.
%
%   SIMZSHIFT(SR,LEN,NP,TP,DELAY,THETAFREQ,UF,UNR) needs 8 input arguments:
%       SR: sampling rate
%       LEN: length of the simulated segments
%       NP: proportion of the added Gaussian noise
%       TP: scaling parameter for the information transfer
%       DELAY: applied delay between unit and LFP
%       THETAFREQ: frequency of LFP oscillation
%       UF: unit frequency
%       UNR: methos for unit generation ('Poisson' or 'tonic')
%
%   See also ZSHIFTRUN5, SIMZSHIFT and SIMZSHIFTRUN.

% Initial parameters
% sr = 1000;      % sampling rate: 1000 Hz
% len = 5;       % segment length: 10 s
% np = 0.1;       % proportion of Gaussian noise to add
% tp = 0.5;       % transfer rate
% delay = 50;    % delay: 50 ms
% unr = 'Poisson';  % method for random unit generation
% thetafreq = 4;  % theta frequency
unitfreq = uf /sr;  % unit frequency

% Random EEG
tim = 0:1/sr:len;
tim = tim(2:end);
eeg = sin(tim*2*pi*thetafreq);
eeg = eeg + rand(size(eeg)) * np;

% Random unit
if isequal(unr,'Poisson')       % random Poisson unit
    lambda = unitfreq;
    r = random('exp',1/lambda,1,10000);     % 1/lambda param. exp.!!!! (MATLAB-ban forditva van...)
    s = cumsum(r);
    pvdisc = unique(ceil(s));
    vdisc = pvdisc(pvdisc<len*sr);   % expected value of length(unit): len*lambda
    if isequal(length(vdisc),length(pvdisc))
        error('Technical error 106')
    end
elseif isequal(unr,'tonic')     % random tonic unit
    vdisc = (1/unitfreq):(1/unitfreq):len*sr;
    vdisc = round(vdisc);
end

% Transfer
su = bhconv(vdisc,len*sr,sr,sr);
d2 = delay / sr * 1000;
if d2 > 0
    sud = [su(d2+1:end) zeros(1,d2)];
else
    sud = [zeros(1,-d2) su(1:end+d2)];
end
eeg = eeg + tp * sud;

% Zshift
[zmaxloc,zmax,hang,hmvl] = zshift(eeg,vdisc,sr);



% -------------------------------------------------------------------------
function zint = bhconv(vdisc,lenu,sr,dsr)
ipunit = zeros(1,lenu);
ipunit(ceil(vdisc/(sr/dsr))) = 1;
wbh = blackmanharris(dsr/5);
wipunit = conv(ipunit,wbh);
ld = length(wipunit) - lenu;
zint = wipunit(floor(ld/2):end-ceil(ld/2)-1);



% ----------------------------------------------------------------------------------
function [zmaxloc,zmax,hang,hmvl] = zshift(eeg,vdisc,sr)

% Filtering EEG
nqf = sr / 2;
flt = fir1(256,[2 8]/nqf,'band');      % bandpass filtering
feeg = filtfilt(flt,1,eeg);

% Hilbert transformation of the eeg
ahee = angle(hilbert(feeg));
aheedeg = ahee * (180 / pi);

% Phase angles - Hilbert
bahee = ahee(vdisc(find(vdisc>0&vdisc<length(eeg))));
n = length(bahee);
ftm = sum(exp(1).^(i*bahee)) / n;    % first trigonometric moment
hang = angle(ftm);   % mean angle
hmvl = abs(ftm);     % mean resultant length

% Shift unit
Z = [];
p = [];
T = [-1*sr:1:1*sr];
for t = T
    vt = vdisc(find(vdisc+t>0&vdisc+t<length(eeg))) + t;  % shift unit
    if isempty(vt)      % Skip segment, if all spikes fall out of range after unit shift
        zmaxloc = NaN;
        zmax = NaN;
        plot(0,0)
        text(0,0,'All spikes out of range after unit shift.','Color','red',...
            'HorizontalAlignment','center');
        return
    end
    bang = ahee(vt);

% Mean resultant length
    n = length(bang);
    ftm = sum(exp(1).^(i*bang)) / n;    % first trigonometric moment
    mrl = abs(ftm);     % mean resultant length
    z = n * (mrl ^ 2);  % Rayleigh's Z statistic
    Z(end+1) = z;
    p(end+1) = exp(1) ^ (-1 * z) * (1 + (2 * z - z ^ 2) / ...
        (4 * n) - (24 * z - 132 * z ^ 2 + 76 * z ^ 3 - 9 * z ^ 4) / (288 * n ^ 2));
end

% Plot result
% plot(T/sr,Z)

sl = 0.005;  % level of significance
siglev = -1 * log(sl);
% hold on
% line(xlim,[siglev siglev],'Color','k','LineStyle','--');

zmax = max(Z);
zmxlc = find(Z==zmax);
% x_lim = xlim;
% y_lim = ylim;
% try
%     text(x_lim(1)+2*(x_lim(2)-x_lim(1))/3,y_lim(1)+2*(y_lim(2)-y_lim(1))/3,num2str(T(zmxlc(1))));
% end
hold off

% Skip file, if phase connection is not significant at the maximum localization
if Z(zmxlc) < siglev
    zmaxloc = NaN;
else
    zmaxloc = T(zmxlc);
end