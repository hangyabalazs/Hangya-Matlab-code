function b_entrysim

% Directories
global DATADIR2
global DATAPATH
inpdir = DATADIR2;
resdir = [DATAPATH 'Entrysim\Data\'];
mm = pwd;
cd(resdir)

% Load
lim1 = 650001;
lim2 = 900000;
ff2 = [inpdir '303_n2_2002_05_06 14_42_15_d.mat'];
load(ff2)
vdisc_orig = vdisc(vdisc>lim1&vdisc<lim2) - lim1;
ff = [inpdir '228_n6_2001_10_01 14_34_24_d.mat'];
load(ff)
eeg_orig = eeg(lim1:lim2);
sr = 10000;

% Create random unit (10 Hz Poisson process)
leneeg = length(eeg_orig);
% frq = 10 / sr;   % in 1/datapoints
% r = random('exp',1/frq,1,10000);     % 1/lambda param. exp.!!!! (MATLAB-ban forditva van...)
% s = cumsum(r);
% psvd = unique(ceil(s));     % 'pseudo vdisc'

% Save
% vdisc = psvd;
vdisc = vdisc_orig;
eeg = eeg_orig;
save random eeg vdisc

% Filter EEG
nqf = sr / 2;
flt = fir1(128,[2/nqf 8/nqf]);
feeg = filtfilt(flt,1,eeg_orig);

% Hilbert transformation of the EEG
ahee = angle(hilbert(feeg));
aheedeg = ahee * (180 / pi);

% Create EEG-dependent unit
zerox = valuecrossing([1:length(ahee)],ahee,0,'up');    % zero-phase corresponds to the theta peaks
psvd_ed = [];
for k = 1:length(zerox)
    lr = logical(round(rand(1)));
%     if lr
        rn = round(rand(1)*5+0.5);
        for kk = 1:rn
            psvd_ed = [psvd_ed zerox(k)+100/1000*sr+kk*20];
        end
%     end
end
% psvd_ed = zerox + 100 / 1000 * sr;      % spikes follow peaks by 100 ms
psvd_ed = round(psvd_ed);
% inx = 1:length(psvd_ed);
% rn = round(rand(size(inx))+0.25);
% rn = ones(size(inx));
% inx2 = inx .* rn;
% inx2 = inx2(inx2>0);
% psvd_ed = psvd_ed(inx2);
figure
plot(feeg)
line([zerox; zerox],[zeros(size(zerox))-1; ones(size(zerox))],'Color','cyan')
line([psvd_ed; psvd_ed],[zeros(size(psvd_ed))-1; ones(size(psvd_ed))],'Color','red')

% Save
vdisc = psvd_ed;
eeg = eeg_orig;
save eeg_dependent eeg vdisc

% Create unit-dependent EEG
gcunit = gaussconv(vdisc_orig,leneeg);
eeg_ud = 0.90 * gcunit + 0.1 * eeg_orig;
figure
subplot(2,1,1)
plot(gcunit)
subplot(2,1,2)
plot(eeg_ud)

% Save
eeg = eeg_ud;
vdisc = vdisc_orig;
save unit_dependent eeg vdisc

cd(mm)



% -------------------------------------------------------------------------
function zint = gaussconv(vdisc,lenu)

ipunit = zeros(1,lenu);
ipunit(vdisc) = 1;
wbh = gausswin(1000);
wipunit = conv(ipunit,wbh);
zint = wipunit(1:lenu);