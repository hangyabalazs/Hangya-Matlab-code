function b_3chconvert
%3CHCONVERT    Converts 3-channel data files to usual data files.
%   3CHCONVERT reads 3-channel data files from source directory and saves
%   two usual type data files to destination directory for each original
%   data file. One contains EEG and MS unit (with M in cell id), the orher
%   contains EEG and HC unit (H in cell id).
%
%   Specify directories via editing the program code!
%
%   See also THRESRUN_3CH.

% Input argument check
error(nargchk(0,0,nargin))

% Directories
global DATADIR
source = [DATADIR '3ch\temp\'];
dest = [DATADIR '3ch_converted\temp\'];
mm = pwd;

% Load & save
cd(source)
files = dir(source);
files = files(3:end);
sf = length(files);
for i = 1:sf
    fname = files(i).name;
    cmps = strread(fname,'%s','delimiter','_');
    fname_MS = [cmps{1} 'M' cmps{2}];
    fname_HC = [cmps{1} 'H' cmps{2}];
    for i = 3:length(cmps)
        fname_MS = [fname_MS '_' cmps{i}];
        fname_HC = [fname_HC '_' cmps{i}];
    end
    
    load(fname)
    HCunit = data(:,1);
    MSunit = data(:,2);
    eeg = data(:,3);
    
    data_old = data;
    data = [eeg MSunit];
    save([dest fname_MS],'data')
    data = [eeg HCunit];
    save([dest fname_HC],'data')
end
cd(mm)