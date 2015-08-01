function [lf lf_lin lf_MI lf_all] = iiocomp
%IIOCOMP   Rate of significant assotiations.
%   [LF LF_LIN LF_MI LF_ALL] = IIOCOMP valculates the number of all 
%   significant linear and non-linear correaltions from all patients
%   (LF_LIN and LF_MI), the number of cases where only linear correaltion
%   was detected (LF) and the number of comparisons (LF_ALL).
%
%   See also IMISIGNIO and ICALLER.

% Read Excel file
global DATAPATH
fn = [DATAPATH 'Ulbert\patients_all.xls'];
headerrows = 0;
[mtx ntx atx] = xlsread(fn);
ntx(1:headerrows,:) = [];
atx(1:headerrows,:) = [];

% Call
lf = zeros(1,26);
lf_lin = zeros(1,26);
lf_MI = zeros(1,26);
lf_all = 0;
for k = 1:26
    disp(k)
    pat = atx{k,1};
    patno = num2str(atx{k,2});
    eg = num2str(atx{k,3});
    nm_rows = atx{k,4};
    nm_cols = atx{k,5};
    esg = [atx{k,6} atx{k,7}];
    
    IM_lin = iio_lin2(pat,patno,eg,nm_rows,nm_cols,esg);
    IM_MI = iio(pat,patno,eg,nm_rows,nm_cols,esg);
    
    lf(k) = length(find(IM_lin&~IM_MI));
    lf_lin(k) = length(find(IM_lin));
    lf_MI(k) = length(find(IM_MI));
    switch k
        case {6,7,8,9,10}
            lf_all = lf_all + numel(IM_lin) - numel(diag(IM_lin)) - size(IM_lin,1) - size(IM_lin,2) + 2;
        case {24,25,26}
            lf_all = lf_all + numel(IM_lin) - numel(diag(IM_lin)) - 6 * size(IM_lin,1) + 9 + 3;
        otherwise
            lf_all = lf_all + numel(IM_lin) - numel(diag(IM_lin));
    end
    
    close all
end

% -------------------------------------------------------------------------
function srIM = iio_lin2(pat,patno,eg,nm_rows,nm_cols,esg)

% Directories
global DATAPATH
inpdir = [DATAPATH 'Ulbert\OITI_' patno '_EEG_' eg '\CCGmap2\'];
inpdir_lin = [DATAPATH 'Ulbert\OITI_' patno '_EEG_' eg '\CCGmap2\'];

% Load MI map
fn = [inpdir 'MIshiftfine.mat'];    % MI map
load(fn)

% Load CCG map
fn = [inpdir_lin 'MIshiftfine.mat'];
load(fn)
sl = [];
fn = [inpdir_lin 'siglev_EEG' eg];
load(fn)
siglev_lin = sl(4,2);   % sig. lev.: 0.0001
rIMlin = rIMax;
rIMlin(rIMax<siglev_lin) = NaN;
rIMlin(isnan(rIMlin)) = 0;

% I/O function 'background'
srIM = squeeze(sum(rIMlin,3));

% -------------------------------------------------------------------------
function srIM = iio(pat,patno,eg,nm_rows,nm_cols,esg)

% Directories
global DATAPATH
inpdir = [DATAPATH 'Ulbert\OITI_' patno '_EEG_' eg '\MImap\'];
inpdir_lin = [DATAPATH 'Ulbert\OITI_' patno '_EEG_' eg '\MImap\'];

% Load MI map
fn = [inpdir 'MIshiftfine.mat'];    % MI map
load(fn)

% Load CCG map
fn = [inpdir_lin 'MIshiftfine.mat'];
load(fn)
sl = [];
fn = [inpdir_lin 'siglev_EEG' eg];
load(fn)
siglev_lin = sl(4,2);   % sig. lev.: 0.0001
rIMlin = rIMax;
rIMlin(rIMax<siglev_lin) = NaN;
rIMlin(isnan(rIMlin)) = 0;

% I/O function 'background'
srIM = squeeze(sum(rIMlin,3));