function somisstim2_call(I)
%SOMISSTIM2_CALL   Calls SOMISSTIM2 for assessment of optical tagging.
%   SOMISSTIM2_CALL calls SOMISSTIM2 for all cells in CellBase (see also
%   CellBase documentation). Output vatiables of SOMISSTIM2 are saved in a
%   result directory.
%
%   SOMISSTIM2_CALL(I) runs only for the index set I.
%
%   See also SOMISSTIM2.

% Pass the control to the user in case of error
dbstop if error

% Directories
global DATAPATH
fs = filesep;
resdir = [DATAPATH 'VIP' fs 'Hindex' fs];

% Load CellBase
loadcb

% Input argument check
nmc = length(CELLIDLIST);
if nargin < 1
    I = 1:nmc;
end

% Call 'somisstim2'
p1 = nan(1,nmc);
p2 = nan(1,nmc);
p3 = nan(1,nmc);
p4 = nan(1,nmc);
I1 = nan(1,nmc);
I2 = nan(1,nmc);
I3 = nan(1,nmc);
I4 = nan(1,nmc);
ff1 = [resdir 'p_val.mat'];
ff2 = [resdir 'I_diff.mat'];
for k = I;
    disp(k)
    try
        [p1(k) I1(k) p2(k) I2(k) p3(k) I3(k) p4(k) I4(k)] = ...
            somisstim3(CELLIDLIST{k});
        save(ff1,'p1','p2','p3','p4')
        save(ff2,'I1','I2','I3','I4')
    catch ME
        disp(['No p values calculated for cell ' num2str(k) '.'])
        disp(ME.message)
    end
end
keyboard