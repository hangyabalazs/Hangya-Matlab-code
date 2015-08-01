function intan_analysis
%INTAN_ANALYSIS   Spike detection.
%   Batch processing for spike detection on Intan data files. Time stamps
%   and waveforms of the detected action potentials are saved.
%
%   See also INTANDISC2.

% Directories
inpdir = 'c:\Balazs\courses\2014_TENSS\data\Brian\test\';
resdir = 'c:\Balazs\courses\2014_TENSS\data\Brian\test\';
dr = dir(inpdir);
dr = dr(3:end);
NumFiles = length(dr);

% Batch
NumTetrodes = 4;
last_ts = 0;
[AllTimeStamps AllWaveForms] = deal(cell(1,NumTetrodes));
for iF = 1:NumFiles
    fn = fullfile(inpdir,dr(iF).name);
    disp([num2str(iF) '/' num2str(NumFiles) '  ' fn])
    [ts1 wf1 t1] = readfa(fn);  % read first half of the file
    [ts2 wf2 t2] = readfb(fn);  % read second half of the file
    t = [t1; t2];   % time
    [ts wf] = deal(cell(1,NumTetrodes));
    for iT = 1:NumTetrodes
        ts{iT} = [ts1{iT} t2(1)+ts2{iT}];  % concatenate time stamps
        wf{iT} = cat(1,wf1{iT},wf2{iT});   % concatenate waveforms
        AllTimeStamps{iT} = [AllTimeStamps{iT} ts{iT}+last_ts];
        AllWaveForms{iT} = cat(1,AllWaveForms{iT},wf{iT});
    end
    last_ts = last_ts + t(end);  % last time stamp
end

% Save
savestr0 = ['save(''' resdir filesep 'TT'];
for iT = 1:NumTetrodes
    TimeStamps = AllTimeStamps{iT}; %#ok<NASGU>
    WaveForms = AllWaveForms{iT}; %#ok<NASGU>
    savestr = [savestr0 num2str(iT) ''',''WaveForms'',''TimeStamps'');'];
    eval(savestr)
end

% -------------------------------------------------------------------------
function [TimeStamps WaveForms t] = readfa(fn)

% Load data
% fn = 'c:\Balazs\courses\2013_TENSS\animal3\recgr13_130613_231048.int';
[t,amps,data,aux] = read_intan_data2a(fn);
data = [data(:,1:3) zeros(size(data,1),1) data(:,4:end)];

% Filter and threshold
thr = 70;
[TimeStamps WaveForms] = intandisc2(data,thr,[]);

% -------------------------------------------------------------------------
function [TimeStamps WaveForms t] = readfb(fn)

% Load data
% fn = 'c:\Balazs\courses\2013_TENSS\animal3\recgr13_130613_231048.int';
[t,amps,data,aux] = read_intan_data2b(fn);
data = [data(:,1:3) zeros(size(data,1),1) data(:,4:end)];

% Filter and threshold
thr = 50;
[TimeStamps WaveForms] = intandisc2(data,thr,[]);