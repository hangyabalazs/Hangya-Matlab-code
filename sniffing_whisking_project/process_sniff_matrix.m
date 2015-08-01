function[fsnm] = process_sniff_matrix(snm,sf,flow)
% fsnm = function process_sniff_matrix(snm,sf,flow)
% independent trials are in columns.
% detrend
% zscore
% filter
% SPR

if nargin<3,
    flow = 100;
end
if nargin < 2,
    sf = 1000;
end
fsnm = nan(size(snm));
for iT = 1:size(snm,2),
    snf_vec = detrend(snm(:,iT));
    snf_vec = zscore(snf_vec);
    fsnm(:,iT) = process_sniffing(snf_vec,sf,flow);
end
