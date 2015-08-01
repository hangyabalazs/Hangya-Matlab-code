function afrerestrictstand_call_layer5
%AFRERESTRICTSTAND_CALL_LAYER5    Calls AFRE_RESTRICT_STAND for a sequence of directories.
%
%   See also AFRE_RESTRICT_STAND_LAYER5.

% Directories
global DATAPATH
% inproot = [DATAPATH 'Hajni_layer_5\disc2\'];
inproot = [DATAPATH 'Hajni\EEGMPO\data\data for EEG_EPSP phase\mat2\'];
dr = dir(inproot);
inpadd = {};
for k = 3:length(dr)
    if dr(k).isdir
        inpadd{end+1} = dr(k).name;
    end
end

% Call AFRE_RESTRICT_STAND
for k = 1:length(inpadd)
    inpdir = [inproot inpadd{k} '\']
%     afre_restrict_stand_layer5(inpdir)
%     L5_spiking_5percent_nocut(inpdir);
    PO_EPSPs_5percent(inpdir);
end