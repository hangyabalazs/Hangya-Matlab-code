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
xlsname = [resdir fs 'tagging1.xls'];   % write results to excel file

% Load CellBase
loadcb

% Input argument check
nmc = length(CELLIDLIST);    %#ok<USENS>
if nargin < 1
    I = 1:nmc;
end

% Call 'Lratio', 'nbisstim', 'nbspikeshapecorr'
feature_names1 = {'Amplitude','Energy'};
feature_names2 = {'WavePC1','Energy'};
for k = I;
    cellid = CELLIDLIST{k};
    disp([num2str(k) '   ' cellid])
    try
        
        % Determine whether there was a stimulation session
        if logical(exist(cellid2fnames(cellid,'StimEvents'),'file'))
            isstim = 1;
        else
            disp(['No ''StimEvents'' file for ' cellid])
            isstim = 0;
        end
        
        % Cluster quality
        [ID_amp Lr_amp] = LRatio(cellid,feature_names1);
        [ID_PC Lr_PC valid_channels] = LRatio(cellid,feature_names2); %#ok<NASGU>
        
        if isstim
            
            % Add 'PulseOn' event if missing
            ST = loadcb(cellid,'STIMSPIKES');
            if isequal(findcellstr(ST.events(:,1),'PulseOn'),0)
                prealignSpikes(CELLIDLIST(k),'FUNdefineEventsEpochs',...
                    @defineEventsEpochs_pulseon,'filetype','stim',...
                    'ifsave',1,'ifappend',1)
            end
            
            % Tagging index
            [Hindex D_KL] = nbisstim(cellid);
            
            % Spike shape correlation
            R = nbspikeshapecorr(cellid);
        
        else
           Hindex = NaN;
           D_KL = NaN;
           R = NaN;
        end
                
        if issave
            
            % Save
            save([resdir 'TAGGING_' regexprep(cellid,'\.','_') '.mat'],...
                'Lr_amp','ID_amp','Lr_PC','ID_PC','valid_channels',...
                'Hindex','D_KL','R')
            
            % Write to Excel
            if isequal(ID_amp,Inf)      % Inf is converted to 65535 in Excel
                ID_amp_xls = {num2str(ID_amp)};
            else
                ID_amp_xls = ID_amp;
            end
            if isequal(ID_PC,Inf)
                ID_PC_xls = {num2str(ID_PC)};
            else
                ID_PC_xls = ID_PC;
            end
            xlswrite(xlsname,CELLIDLIST(k),'sheet1',['A' num2str(k)])
            xlswrite(xlsname,Lr_amp,'sheet1',['B' num2str(k)])
            xlswrite(xlsname,ID_amp_xls,'sheet1',['C' num2str(k)])
            xlswrite(xlsname,Lr_PC,'sheet1',['D' num2str(k)])
            xlswrite(xlsname,ID_PC_xls,'sheet1',['E' num2str(k)])
            xlswrite(xlsname,Hindex,'sheet1',['F' num2str(k)])
            xlswrite(xlsname,D_KL,'sheet1',['G' num2str(k)])
            xlswrite(xlsname,R,'sheet1',['H' num2str(k)])
        end
        
    catch ME
        disp(['Something went wrong for cell ' num2str(k) '.'])
        disp(ME.message)
    end
end
% keyboard