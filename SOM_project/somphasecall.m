function somphasecall
%SOMPHASECALL   Calls SOMPHASE.
%   SOMPHASECALL calls SOMPHASE for a sequence of data files. See SOMPHASE
%   for more details.
%
%   See also SOMPHASE.

% Directories
global DATAPATH
inpdir = '\\science\Kepecs\recordings\for Balazs\new\';
resdir = [DATAPATH 'SOM\Phase\'];
mm = pwd;
cd(resdir)

% Phase histogram
dr = dir(inpdir);
for k = 3:length(dr)
    if dr(k).isdir
        
        % Load EEG
        inpadd = dr(k).name;
        drr = [inpdir inpadd];
        [fl0 fl] = filelist(drr);
        fst = find(strncmp(fl,'CSC',3));
        eegname = fl{fst(1)};
        eegnm = eegname(1:end-4);
        load([inpdir inpadd '\' eegname])
        eval(['frs = ' eegnm '_SampleFrequencies;'])    % sampling rate
        if max(frs) > 1020
            continue
        end
        eval(['eeg = ' eegnm '_Samples;'])  % EEG
        eeg = eeg(:);
        
        eval(['tseeg = ' eegnm '_TimeStamps;'])    % EEG time stamps
        sr = 1/(tseeg(2)-tseeg(1))*1000000*512;
        dt = 1 / sr;
        time = (0:length(eeg)-1)*dt+CSC1_TimeStamps(1)/1000000;
        
        % Load unit
        fst = find(strncmp(fl,'TT',2));
        for cells = 1:length(fst)
            unitname = fl{fst(cells)};
            unitnm = unitname(1:end-4);
            load([inpdir inpadd '\' unitname])
            
            % Phase calculation
            [hang, hmvl, ftm, bahee] = somphase(eeg,TS,sr,time,3,10);   % theta phase
            cd('theta')
            fnm = ['THETA_PHASE_' inpadd '_' eegnm '_' unitnm '.fig'];
            saveas(gcf,fnm)
            fnm = ['THETA_PHASE_' inpadd '_' eegnm '_' unitnm '.mat'];
            save(fnm,'hang','hmvl','ftm','bahee')
            cd ..
            
            [hang, hmvl, ftm, bahee] = somphase(eeg,TS,sr,time,20,35);  % delta phase
            cd('beta')
            fnm = ['BETA_PHASE_' inpadd '_' eegnm '_' unitnm '.fig'];
            saveas(gcf,fnm)
            fnm = ['BETA_PHASE_' inpadd '_' eegnm '_' unitnm '.mat'];
            save(fnm,'hang','hmvl','ftm','bahee')
            cd ..
        end
    end
end
    
    
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
    