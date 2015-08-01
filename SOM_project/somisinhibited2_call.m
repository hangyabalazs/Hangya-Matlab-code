function somisinhibited2_call(I)
%SOMISINHIBITED2_CALL   Calls SOMISINHIBITED2 for testing inhibition.
%   SOMISINHIBITED2_CALL calls SOMISINHIBITED2 for all cells in CellBase
%   (see also CellBase documentation). Output variables of SOMISINHIBITED2
%   are saved in a result directory.
%
%   SOMISINHIBITED2_CALL(I) runs only for the index set I.
%
%   See also SOMISINHIBITED2.

% Pass the control to the user in case of error
dbstop if error

% Directories
global DATAPATH
fs = filesep;
resdir = [DATAPATH 'SOM' fs 'somisinhibited2' fs];

% Load CellBase
loadcb

% Input argument check
nmc = length(CELLIDLIST);
if nargin < 1
    I = 1:nmc;
end

% Call 'somisinhibited'
p_val = nan(1,nmc);
inhibition_start = nan(1,nmc);
inhibition_end = nan(1,nmc);
inhibition_time = nan(1,nmc);
ff1 = [resdir 'p_val.mat'];
for k = I;
    disp(k)
    cellid = CELLIDLIST{k};
    try
        [p_val(k) inhibition_start(k) inhibition_end(k) inhibition_time(k) H] ...
            = somisinhibited2(cellid);
%         save(ff1,'p_val','inhibition_start','inhibition_end','inhibition_time')   % save
        ff2 = [resdir 'SDF_' cellid '.fig'];
%         saveas(H,ff2)
        close all
    catch ME
        disp(['No p value calculated for cell ' num2str(k) '.'])
        disp(ME.message)
    end
end
keyboard