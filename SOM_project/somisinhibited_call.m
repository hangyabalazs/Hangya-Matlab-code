function somisinhibited_call(I)
%SOMISINHIBITED_CALL   Calls SOMISINHIBITED for testing inhibition.
%   SOMISINHIBITED_CALL calls SOMISINHIBITED for all cells in CellBase (see
%   also CellBase documentation). Output vatiables of SOMISINHIBITED are
%   saved in a result directory.
%
%   SOMISINHIBITED_CALL(I) runs only for the index set I.
%
%   See also SOMISINHIBITED.

% Pass the control to the user in case of error
dbstop if error

% Directories
global DATAPATH
fs = filesep;
resdir = [DATAPATH 'SOM' fs 'somisinhibited' fs];

% Load CellBase
loadcb

% Input argument check
nmc = length(CELLIDLIST);
if nargin < 1
    I = 1:nmc;
end

% Call 'somisinhibited'
p_val = nan(1,nmc);
ff1 = [resdir 'p_val.mat'];
for k = I;
    disp(k)
    try
        p_val(k) = somisinhibited(CELLIDLIST{k});
        save(ff1,'p_val')
    catch ME
        disp(['No p value calculated for cell ' num2str(k) '.'])
        disp(ME.message)
    end
end
keyboard