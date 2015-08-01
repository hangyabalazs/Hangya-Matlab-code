function [rIMaxLoc, rIMax] = imishiftcall_N2(data,sr,delay)
%IMISHIFTCALL_N2   Calls IMISHIFT.
%   IMISHIFTCALL_N2 is a caller function for IMISHIFT calculating time-shifted
%   mutual information.
%   It plots and saves IMISHIFT output.
%
%   [RIML RIM] = IMISHIFTCALL_N2(DATA,SR,D) calls IMISHIFT with input
%   parameters DATA delayed by D, and SR (sampling rate). It returns
%   IMISHIFT output parameters: RIM (maximal mutual information) and RIML
%   (mutual information maximum location).
%
%   Note: prefiltering data (40 Hz, lowpass) is recommended (see IFILTER)!
%
%   See also IMISHIFT and IFILTER.

% Directories
global DATAPATH
eg = '90c';
resdir = [DATAPATH 'Ulbert\OITI_N2_EEG_' eg '\MImap\'];
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
        if x ~= 5 && x ~= 17 && x ~= 18 && y ~= 5 && y ~= 17 && y ~= 18
            [rIMaxLoc(x,y,:) rIMax(x,y,:)] = imishift(ch1,ch2,sr);  % 5, 26 and 27 channels are bad
            [rIMaxLoc(y,x,:) rIMax(y,x,:)] = imishift(ch2,ch1,sr);
        end
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