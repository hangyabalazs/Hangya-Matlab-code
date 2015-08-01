function [rIMaxLoc, rIMax] = imishiftcall_io(data,sr,delay)
%IMISHIFTCALL_IO   Calls IMISHIFT_IO.
%   IMISHIFTCALL_IO is a caller function for IMISHIFT_IO calculating
%   time-shifted mutual information. It plots and saves IMISHIFT_IO output.
%
%   [RIML RIM] = IMISHIFTCALL_IO(DATA,SR,D) calls IMISHIFT_IO with input
%   parameters DATA delayed by D, and SR (sampling rate). It returns
%   IMISHIFT_IO output parameters: RIM (maximal mutual information) and
%   RIML (mutual information maximum location).
%
%   Note: prefiltering data (0.1-40 Hz) is recommended (see IFILTER3)!
%
%   See also IMISHIFT_IO and IFILTER3.

% Directories
global DATAPATH
eg = '90a';
resdir = [DATAPATH 'Ulbert\OITI_N2_EEG_' eg '\MImap\'];
resdir_io = [DATAPATH 'Ulbert\OITI_N2_EEG_' eg '\IO\'];
mm = pwd;
cd(resdir)  % check whether directories exist
cd(resdir_io)

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
        [rIMaxLoc(x,y,:) rIMax(x,y,:) Hsd Hse fcn fsd fsem fy] = imishift_io(ch1,ch2,sr);
        cd(resdir_io)
        fnm = ['IOsd_' num2str(delay) '_' num2str(x) '_' num2str(y) '.fig'];
        saveas(Hse,fnm)
        fnm = ['IOse_' num2str(delay) '_' num2str(x) '_' num2str(y) '.fig'];
        saveas(Hse,fnm)
        fnm = ['IO_' num2str(delay) '_' num2str(x) '_' num2str(y) '.mat'];
        save(fnm,'fcn','fsd','fsem','fy')
        
        [rIMaxLoc(y,x,:) rIMax(y,x,:) Hsd Hse fcn fsd fsem fy] = imishift_io(ch2,ch1,sr);
        cd(resdir_io)
        fnm = ['IOsd_' num2str(delay) '_' num2str(y) '_' num2str(x) '.fig'];
        saveas(Hse,fnm)
        fnm = ['IOse_' num2str(delay) '_' num2str(y) '_' num2str(x) '.fig'];
        saveas(Hse,fnm)
        fnm = ['IO_' num2str(delay) '_' num2str(y) '_' num2str(x) '.mat'];
        save(fnm,'fcn','fsd','fsem','fy')
        close all
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
% cd(resdir)
% save(['MIshiftmaps_EEG' eg '_' num2str(delay)],'rIMaxLoc','rIMax')
% cd(mm)