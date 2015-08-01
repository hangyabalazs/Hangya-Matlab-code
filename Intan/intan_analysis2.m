function intan_analysis2
%INTAN_ANALYSIS   Spike detection.
%   Batch processing for spike detection on Intan data files. Time stamps
%   and waveforms of the detected action potentials are saved.
%
%   See also INTANDISC2.

% Directories
inpdir = 'c:\Balazs\courses\2013_TENSS\animal3\balazs_test\';
resdir = 'C:\Balazs\courses\2013_TENSS\animal3\balazs_test_processed\';
dr = dir(inpdir);
dr = dr(3:end);
NumFiles = length(dr);

% Batch
NumTetrodes = 4;
last_ts = 0;
[AllTimeStamps AllWaveForms] = deal(cell(1,NumTetrodes));
TTLs = [];
for iF = 1:NumFiles
    fn = fullfile(inpdir,dr(iF).name);
    disp([num2str(iF) '/' num2str(NumFiles) '  ' fn])
    [ts1 wf1 ttls1 t1] = readfa(fn,70);  % read first half of the file
    [ts2 wf2 ttls2 t2] = readfb(fn,70);  % read second half of the file
    t = [t1; t2];   % time
    ttls = [ttls1; ttls2];   % time
    [ts wf] = deal(cell(1,NumTetrodes));
    for iT = 1:NumTetrodes
        ts{iT} = [ts1{iT} t2(1)+ts2{iT}];  % concatenate time stamps
        wf{iT} = cat(1,wf1{iT},wf2{iT});   % concatenate waveforms
        AllTimeStamps{iT} = [AllTimeStamps{iT} ts{iT}+last_ts];
        AllWaveForms{iT} = cat(1,AllWaveForms{iT},wf{iT});
    end
    TTLs = [TTLs; ttls+last_ts];  %#ok<AGROW>   % concatenate TTLs
    last_ts = last_ts + t(end);  % last time stamp
end

% Save
savestr0 = ['save(''' resdir filesep 'TT'];
for iT = 1:NumTetrodes
    TimeStamps = AllTimeStamps{iT}; %#ok<NASGU>
    WaveForms = AllWaveForms{iT}; %#ok<NASGU>
    savestr = [savestr0 num2str(iT) ''',''WaveForms'',''TimeStamps'');'];
    eval(savestr)
    save([resdir filesep 'position_timestamps.mat'],'TTLs')
end

% -------------------------------------------------------------------------
function [TimeStamps WaveForms TTLs t] = readfa(fn,thr)

% Load data
% fn = 'c:\Balazs\courses\2013_TENSS\animal3\recgr13_130613_231048.int';
[t,amps,data,aux] = read_intan_data2a(fn);

% Filter and threshold
[TimeStamps WaveForms] = intandisc2(data,thr,[]);

% TTLs
TTLs = t(diff(double(aux(:,6)))~=0);

% -------------------------------------------------------------------------
function [TimeStamps WaveForms TTLs t] = readfb(fn,thr)

% Load data
% fn = 'c:\Balazs\courses\2013_TENSS\animal3\recgr13_130613_231048.int';
[t,amps,data,aux] = read_intan_data2b(fn);

% Filter and threshold
[TimeStamps WaveForms] = intandisc2(data,thr,[]);

% TTLs
TTLs = t(diff(double(aux(:,6)))~=0);