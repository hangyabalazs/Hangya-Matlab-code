function vipisinfluenced_call(I)
%VIPISINFLUENCED_CALL   Calls VIPISINFLUENCED2 for testing inhibition.
%   VIPISINFLUENCED_CALL calls VIPISINFLUENCED2 for all cells in CellBase
%   (see also CellBase documentation). Output variables of VIPISINFLUENCED2
%   are saved in a result directory.
%
%   VIPISINFLUENCED_CALL(I) runs only for the index set I.
%
%   See also VIPISINFLUENCED2.

% Pass the control to the user in case of error
dbstop if error

% Directories
global DATAPATH
fs = filesep;
resdir = [DATAPATH 'VIP' fs 'vipisinfluenced11_gonogo2_p' fs];

% Load CellBase
loadcb

% Input argument check
nmc = length(CELLIDLIST);
if nargin < 1
    I = 1:nmc;
end

% Call 'vipisinfluenced'
p_inh = nan(nmc,1);
p_act = nan(nmc,1);
baseline = nan(nmc,1);
minvalue = nan(nmc,1);
inhibition_start = nan(nmc,1);
inhibition_end = nan(nmc,1);
inhibition_peak = nan(nmc,1);
inhibition_time = nan(nmc,1);
maxvalue = nan(nmc,1);
activation_start = nan(nmc,1);
activation_end = nan(nmc,1);
activation_peak = nan(nmc,1);
activation_time = nan(nmc,1);
ff1 = [resdir 'p_val.mat'];
for k = I;
    disp(k)
    cellid = CELLIDLIST{k};
    try
%         [p_inh(k) inhibition_start(k) inhibition_end(k) inhibition_time(k) ...
%             p_act(k) activation_start(k) activation_end(k) activation_time(k) H Hr] ...
%             = vipisinfluenced2(cellid);
%         [baseline(k) p_inh(k) minvalue(k) inhibition_start(k) inhibition_end(k) inhibition_peak(k) inhibition_time(k) ...
%             p_act(k) maxvalue(k) activation_start(k) activation_end(k) activation_peak(k) activation_time(k) H Hr] = ...
%             vipisinfluenced3_burston(cellid);
%         [baseline(k) p_inh(k) minvalue(k) inhibition_start(k) inhibition_end(k) inhibition_peak(k) inhibition_time(k) ...
%             p_act(k) maxvalue(k) activation_start(k) activation_end(k) activation_peak(k) activation_time(k) H Hr] = ...
%             vipisinfluenced3(cellid);
        [baseline(k) p_inh(k) minvalue(k) inhibition_start(k) inhibition_end(k) inhibition_peak(k) inhibition_time(k) ...
            p_act(k) maxvalue(k) activation_start(k) activation_end(k) activation_peak(k) activation_time(k) H Hr] = ...
            vipisinfluenced3c(cellid);
%         [baseline(k) p_inh(k) minvalue(k) inhibition_start(k) inhibition_end(k) inhibition_peak(k) inhibition_time(k) ...
%             p_act(k) maxvalue(k) activation_start(k) activation_end(k) activation_peak(k) activation_time(k) H Hr] = ...
%             vipisinfluenced3b_A1(cellid);
%         save(ff1,'p_inh','p_act','inhibition_start','inhibition_end','inhibition_time',...
%             'p_act','activation_start','activation_end','activation_time')   % save
        save(ff1,'baseline','minvalue','p_inh','p_act','inhibition_start','inhibition_end','inhibition_peak','inhibition_time',...
            'maxvalue','p_act','activation_start','activation_end','activation_peak','activation_time')   % save
        ff2 = [resdir 'SDF_' cellid '.fig'];
        saveas(H,ff2)
        ff3 = [resdir 'RASTER_' cellid '.fig'];
        saveas(Hr,ff3)
        close all
    catch ME
        disp(['No p value calculated for cell ' num2str(k) '.'])
        disp(ME.message)
    end
end
keyboard