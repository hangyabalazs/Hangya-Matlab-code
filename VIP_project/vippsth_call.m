function vippsth_call(I)
%VIPPSTH_CALL   Calls VIPPSTH.
%   VIPPSTH_CALL calls VIPPSTH for all cells in CellBase (see also CellBase
%   documentation). 
%
%   VIPPSTH_CALL(I) runs only for the index set I.
%
%   See also VIPPSTH.

% Pass the control to the user in case of error
dbstop if error

% Directories
global DATAPATH
fs = filesep;
resdir = [DATAPATH 'VIP' fs 'vippsth2_A1' fs];

% Load CellBase
loadcb

% Input argument check
nmc = length(CELLIDLIST);
if nargin < 1
    I = 1:nmc;
end

% Call 'vippsth'
for k = I;
    disp(k)
    cellid = CELLIDLIST{k};
    try
        H = vippsth_A1(cellid);
        ff2 = [resdir 'PSTH_' cellid '.fig'];
        saveas(H,ff2)
        close(H)
    catch ME
        disp(['No PSTH calculated for cell ' num2str(k) '.'])
        disp(ME.message)
    end
end
keyboard