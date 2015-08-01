function [rIMaxLoc, rIMax] = imishiftcall_lowres(data,sr,delay)
%IMISHIFTCALL_LOWRES   Calls IMISHIFT_LOWRES.
%   IMISHIFTCALL_LOWRES is a caller function for IMISHIFT_LOWRES
%   calculating time-shifted mutual information. It plots and saves
%   IMISHIFT_LOWRES output.
%
%   [RIML RIM] = IMISHIFTCALL_LOWRES(DATA,SR,D) calls IMISHIFT_LOWRES with
%   input parameters DATA delayed by D, and SR (sampling rate). It returns
%   IMISHIFT output parameters: RIM (maximal mutual information) and RIML
%   (mutual information maximum location).
%
%   Note: prefiltering data (0.1-40 Hz, bandpass) is recommended (see
%   IFILTER3)!
%
%   See also IMISHIFT and IFILTER.

% Directories
global DATAPATH
eg = '177';
resdir = [DATAPATH 'Ulbert\OITI_31_EEG_' eg '_long\MImap\'];
mm = pwd;
cd(resdir)

% Main
chno = size(data,2);    % prefiltering data (40 Hz, lowpass) is recommended (see IFILTER)!
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
        [rIMaxLoc(x,y,:) rIMax(x,y,:)] = imishift_lowres(ch1,ch2,sr);
        [rIMaxLoc(y,x,:) rIMax(y,x,:)] = imishift_lowres(ch2,ch1,sr);
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
save(['MIshiftmaps_EEG' eg 'long_' num2str(delay)],'rIMaxLoc','rIMax')
cd(mm)