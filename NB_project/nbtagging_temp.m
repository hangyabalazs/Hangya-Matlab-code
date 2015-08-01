function nbtagging(I,issave)
%NBTAGGING   Cluster quality and tagging index.
%   NBTAGGING calculates cluster quality measures, tagging index and spike
%   shape correlation between spontaneous and light-evoked action
%   potentials for all cells in CellBase (see also CellBase documentation).
%   It calls LRATIO, NBTAGGING and NBSPIKESHAPECORR.
%
%   NBTAGGING(I) runs only for the index set I.
%
%   NBTAGGING(I,IS) accepts a second, logical input argument determining
%   whether to save the results of the analysis. Default beahvior is not
%   saving.
%
%   See also LRATIO, NBTAGGING and NBSPIKESHAPECORR.

% Pass the control to the user in case of error
dbstop if error

% Input argument check
error(nargchk(1,2,nargin))
if nargin < 2
    issave = false;
end

% Directories
global DATAPATH
fs = filesep;
resdir = [DATAPATH 'NB' fs 'tagging' fs];
% dsr = datestr(now);
% dsr = regexprep(dsr,':','_');
% ff1 = [resdir 'cluster_quality_' dsr '.mat'];
% ff2 = [resdir 'Hindex_' dsr '.mat'];
% ff3 = [resdir 'spikeshape_correlation_' dsr '.mat'];
xlsname = [resdir fs 'tagging.xls'];   % write results to excel file

% Load CellBase
loadcb

% Input argument check
nmc = length(CELLIDLIST);    %#ok<USENS>
if nargin < 1
    I = 1:nmc;
end

% Hindex for pulseon, frequency-wise; for burston (effect of first stim);
% L-ratio, ID; spike shape correlation; latency, jitter; reliability
% (frequency-wise also); distance from light-evoked noise spikes

% first round: ID, L-ratio, Hindex for pulseon (downsample if two many), waveform correlation

% Call
feature_names1 = {'Amplitude','Energy'};
feature_names2 = {'WavePC1','Energy'};
Lr_amp = nan(nmc,1);
ID_amp = nan(nmc,1);
Lr_PC = nan(nmc,1);
ID_PC = nan(nmc,1);
valid_channels = nan(nmc,4);
Hindex = nan(nmc,1);
D_KL = nan(nmc,1);
R = nan(nmc,1);
for k = I;
    disp(k)
    cellid = CELLIDLIST{k};
    try
        % Add 'PulseOn' event if missing
        ST = loadcb(cellid,'STIMSPIKES');
        if isequal(findcellstr(ST.events(:,1),'PulseOn'),0)
            prealignSpikes(CELLIDLIST(k),'FUNdefineEventsEpochs',...
                @defineEventsEpochs_pulseon,'filetype','stim',...
                'ifsave',1,'ifappend',1)
        end
        
        % Cluster quality
        [ID_amp(k) Lr_amp(k)] = LRatio(cellid,feature_names1);
        [ID_PC(k) Lr_PC(k) valid_channels(k,:)] = LRatio(cellid,feature_names2);
                
        % Tagging index
        [Hindex(k) D_KL(k)] = nbisstim(cellid);
                
        % Spike shape correlation
        R(k) = nbspikeshapecorr(cellid);
                
        if issave
            
            % Save
%             save(ff1,'Lr_amp','ID_amp','Lr_PC','ID_PC')
%             save(ff2,'Hindex','D_KL')
%             save(ff3,'R')
            save([resdir 'TAGGING_' regexprep(cellid,'\.','_') '.mat'],...
                'Lr_amp','ID_amp','Lr_PC','ID_PC','Hindex','D_KL','R')
            
            % Write to Excel
            xlswrite(xlsname,CELLIDLIST(k),'sheet1',['A' num2str(k)])
            xlswrite(xlsname,Lr_amp(k),'sheet1',['B' num2str(k)])
            xlswrite(xlsname,ID_amp(k),'sheet1',['C' num2str(k)])
            xlswrite(xlsname,Lr_PC(k),'sheet1',['D' num2str(k)])
            xlswrite(xlsname,ID_PC(k),'sheet1',['E' num2str(k)])
            xlswrite(xlsname,Hindex(k),'sheet1',['F' num2str(k)])
            xlswrite(xlsname,D_KL(k),'sheet1',['G' num2str(k)])
            xlswrite(xlsname,R(k),'sheet1',['H' num2str(k)])
        end
        
    catch ME
        disp(['Something went wrong for cell ' num2str(k) '.'])
        disp(ME.message)
    end
end
% keyboard