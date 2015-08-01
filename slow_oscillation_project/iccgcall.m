function [rIMaxLoc, rIMax] = iccgcall(data,sr,delay)
%ICCGCALL   Calls ICCG.
%   ICCGCALL is a caller function for ICCG calculating linear
%   cross-correlation. It plots and saves ICCG output.
%
%   [RIML RIM] = ICCGCALL(DATA,SR,D) calls ICCG with input parameters DATA
%   delayed by D, and SR (sampling rate). It returns ICCG output
%   parameters: RIM (maximal correlation) and RIML (maximum location).
%
%   Note: prefiltering data (0.1-40 Hz) is recommended (see IFILTER3)!
%
%   See also ICCG and IFILTER3.

% Directories
global DATAPATH
eg = '90a';
resdir = [DATAPATH 'Ulbert\OITI_N2_EEG_' eg '\CCGmap2\'];
mm = pwd;
try
    cd(resdir)
catch
    mkdir(resdir)
    cd(resdir)
end

% Filter
data = lfilter(data,1000);

% Main
chno = size(data,2);    % prefiltering data (0.1-40 Hz) is recommended (see IFILTER3)!
lene = floor((size(data,1)-delay)/sr);
rIMaxLoc = zeros(chno,chno,lene);
rIMax = zeros(chno,chno,lene);
delay = max(1,delay);
for x = 1:chno
    for y = x+1:chno
        disp([num2str(x) ' ' num2str(y)])
        ch1 = data(delay:end,x)';
        ch2 = data(delay:end,y)';
%         profile off
%         profile on -detail builtin -timer real
%         tic
        [rIMaxLoc(x,y,:) rIMax(x,y,:)] = iccg(ch1,ch2,sr);
        [rIMaxLoc(y,x,:) rIMax(y,x,:)] = iccg(ch2,ch1,sr);
%         toc
%         profile viewer
    end
end

% Plot
figure;
imagesc(nanmean(rIMaxLoc,3))
figure;
imagesc(nanmean(rIMax,3))

% Save
cd(resdir)
save(['MIshiftmaps_EEG' eg '_' num2str(delay)],'rIMaxLoc','rIMax')
cd(mm)

% -------------------------------------------------------------------------
function data = lfilter(data,sr)
%IFILTER3   Filters multichannel EEG.
%   FD = IFILTER3(DATA,SR) filters DATA sampled on SR between 0.5 and 4 Hz
%   using a 4096 order lowpass FIR filter. Filtered data is returned in FD.
%
%   See also FIR1 and FILTFILT.

% Construct filter
nqf = sr / 2;      % Nyquist frequency
flt = fir1(2*4096,[0.5 4]/nqf,'band');      % bandpass filtering between 0.1 and 40 Hz

% Filter EEG
chno = size(data,2);
feeg = zeros(size(data));
for x = 1:chno
    feeg(:,x) = filtfilt(flt,1,data(:,x));
end

% Return filtered data
data = feeg;