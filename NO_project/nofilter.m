function nofilter
%NOFILTER   Filtering EEG.
%   NOFILTER applies zero phase shift FIR filtering on EEG files of the NO
%   project.
%
%   See also FIR1 and FILTFILT.

% Directories
global DATADIR
global DATAPATH
inpdir = [DATADIR 'NO_matfiles\'];
resdir = [DATAPATH 'NO\Wavelet\'];
mm = pwd;
cd(resdir)

% Filelist
[files files_short] = filelist(inpdir);
sf = length(files_short);

% Progress indicator
wb = waitbar(0,'Running NOSPECT...','Position',[360 250 275 50]);
global WB
WB(end+1) = wb;

% Main
sr = 5000;     % sampling rate
nqf = sr / 2;   % Nyquist frequency
for o = 2:sf
    fname = files(o).name;
    cmps = strread(fname(1:end-4),'%s','delimiter','_');
    titlestr = [cmps{1} ' ' cmps{2}];
    
    load([inpdir fname])
    eeg = hEEG.values;
    eeg2 = eeg(1:5:end)';
    flt = fir1(4096,2/500,'low');
    feeg = filtfilt(flt,1,eeg2)';
    flt = fir1(4096,[2 8]/500,'band');
    feeg_theta = filtfilt(flt,1,eeg2)';
    figure
    plot(eeg2)
    hold on
    plot(feeg,'r')
    
    
    
%     fn = [fname(1:end-4) '_EEGWAVELET.jpg'];
%     saveas(gcf,fn);
    close all
    waitbar(o/sf)
end
close(wb)
cd(mm)

% -------------------------------------------------------------------------
function [files2 files2_short] = filelist(inpdir)

% List of filenames
files = dir(inpdir);
files = files(3:end);
files2 = struct('name',[],'date',[],'bytes',[],'isdir',[],'datenum',[]);
files2_short = {};
for i = 1:length(files)
    if ~files(i).isdir
        files2(end+1) = files(i);
        files2_short{end+1} = files(i).name;
    end
end
files2 = files2(2:end);