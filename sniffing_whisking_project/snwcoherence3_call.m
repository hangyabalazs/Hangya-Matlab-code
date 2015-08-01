function snwcoherence3_call
%SNWCOHERENCE3_CALL   Caller function for SNWCOHERENCE3.
%   SNWCOHERENCE3_CALL calls SNWCOHERENCE3 to perform phase analysis on
%   sniffing-whisking data.
%
%   See also SNWCOHERENCE3.

% Import session data from Excel
global DATAPATH
tblfile = [DATAPATH 'SniffWhisk\data_selection_coherence.xlsx'];
[tbl0 tbl] = xlsread(tblfile);

% Call 'snwdisc'
numSessions = size(tbl,1);
for iS = 1:numSessions
    rrat = tbl{iS,1};
    rdate = tbl{iS,2};
    csc1 = tbl{iS,3};
    csc2 = tbl{iS,4};
    msgn = tbl0(iS,1);
    wl = tbl0(iS,3);
    epoch_start = tbl0(iS,4);
    epoch_end = tbl0(iS,5);
    
    snwcoherence3(rrat,rdate,csc1,csc2,msgn,wl,epoch_start,epoch_end)
    close all
end