function [TIM SNR] = isnrdist_call
%ISNRDIST_CALL   Caller function for ISNRDIST.
%
%   See also ISNRDIST.

% Read Excel file
global DATAPATH
fn = [DATAPATH 'Ulbert\patients_all.xls'];
headerrows = 0;
[mtx ntx atx] = xlsread(fn);
ntx(1:headerrows,:) = [];
atx(1:headerrows,:) = [];

% Call
TIM = [];
SNR = [];
for k = 24:24
%     k
    pat = atx{k,1};
    patno = num2str(atx{k,2});
    eg = num2str(atx{k,3});
    
    [pTIM pSNR] = isnrdist(pat,patno,eg);
    TIM = [TIM; pTIM];
    SNR = [SNR; pSNR];
    close all
end