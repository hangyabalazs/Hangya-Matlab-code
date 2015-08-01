function nbspikeshapecorr_call(I,issave)
%NBSPIKESHAPECORR_CALL   Call SPIKESHAPECORR.
%   NBSPIKESHAPECORR_CALL calculates  spike shape correlation between
%   spontaneous and light-evoked action potentials for all cells in
%   CellBase (see also CellBase documentation). It calls SPIKESHAPECORR.
%
%   SPIKESHAPECORR_CALL(I) runs only for the index set I.
%
%   SPIKESHAPECORR_CALL(I,IS) accepts a second, logical input argument
%   determining whether to save the results of the analysis.
%
%   See also SPIKESHAPECORR and TAGGING.

% Pass the control to the user in case of error
dbstop if error

% Input argument check
error(nargchk(0,2,nargin))
if nargin < 2
    issave = true;
end

% Directories
global DATAPATH
fs = filesep;
resdir = [DATAPATH 'NB' fs 'tagging2' fs];
xlsname = [resdir fs 'spikeshapecorr.xls'];   % write results to excel file

% Load CellBase
loadcb

% Input argument check
nmc = length(CELLIDLIST);    %#ok<USENS>
if nargin < 1
    I = 1:nmc;
end

% Call 'spikeshapecorr'
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
        
        if isstim
            
            % Spike shape correlation
            R = spikeshapecorr(cellid);
        
        else
           R = NaN;
        end
                
        if issave
            
            % Save
            save([resdir 'SPIKESHAPECORR_' regexprep(cellid,'\.','_') '.mat'],'R')
            
            % Write to Excel
            xlswrite(xlsname,CELLIDLIST(k),'sheet1',['A' num2str(k)])
            xlswrite(xlsname,R,'sheet1',['B' num2str(k)])
        end
        
    catch ME
        disp(['Something went wrong for cell ' num2str(k) '.'])
        disp(ME.message)
    end
end
% keyboard