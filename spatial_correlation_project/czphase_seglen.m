function czphase_seglen
%CZPHASE_SEGLEN   Length of recording segment.
%   CZPHASE_SEGLEN displays the lengths of recording segments for place
%   cell-interneuron pair recordings in the Command Window.
%
%   See also CZRIPPLE.

% Directories
global DATAPATH
inpdir_eeg = [DATAPATH 'Czurko\czurko_EEG\'];

% Load
xlsname = [inpdir_eeg 'EEG3.xls'];
[ntz mtz] = xlsread(xlsname,'all_unique');
sf = size(mtz,1);   % number of pairs
for o = 1:sf
    
    fn = [inpdir_eeg 'EEG_' mtz{o,3} '_' mtz{o,1} '.mat'];
    load(fn)        % load EEG
    eval(['Eeg = ' mtz{o,3} ';']);
    eval(['clear ' mtz{o,3}]);
    eeg = Eeg.values;
    sr = 1 / Eeg.interval;
    eeg_start = Eeg.start;
%     eeg_start = 0;
    eeg_end = eeg_start + (Eeg.length - 1) * Eeg.interval;
    eeg_times = eeg_start:Eeg.interval:eeg_end;
    eeg_end / 60
    
end