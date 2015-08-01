function snwdisc_call
%SNWDISC_CALL   Caller function for SNWDISC.
%   SNWDISC_CALL calls SNWDISC to perform threshold based discimination on
%   whisking data.
%
%   See also SNWDISC.

% Import session data from Excel
global DATAPATH
tblfile = [DATAPATH 'SniffWhisk\data_selection.xlsx'];
[tbl0 tbl] = xlsread(tblfile);

% Call 'snwdisc'
numSessions = size(tbl,1);
for iS = 7:numSessions
    rrat = tbl{iS,1};
    rdate = tbl{iS,2};
    csc1 = tbl{iS,3};
    csc2 = tbl{iS,4};
    msgn = tbl0(iS,1);
    snwdisc(rrat,rdate,csc1,csc2,msgn)
end