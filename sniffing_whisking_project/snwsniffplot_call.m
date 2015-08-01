function snwsniffplot_call
%SNWSNIFFPLOT_CALL   Caller function for SNWSNIFFPLOT.
%   SNWSNIFFPLOT_CALL calls SNWSNIFFPLOT to perform phase analysis on
%   sniffing-whisking data.
%
%   See also SNWSNIFFPLOT.

% Import session data from Excel
global DATAPATH
tblfile = [DATAPATH 'SniffWhisk\data_selection.xlsx'];
[tbl0 tbl] = xlsread(tblfile);

% Call 'snwdisc'
numSessions = size(tbl,1);
for iS = 1:numSessions
    rrat = tbl{iS,1};
    rdate = tbl{iS,2};
    csc1 = tbl{iS,3};
    csc2 = tbl{iS,4};
    msgn = tbl0(iS,1);
    thr = tbl0(iS,2);
    wl = tbl0(iS,3);
    snwsniffplot2(rrat,rdate,csc1,csc2,msgn,thr,wl)
%     snwcoherence2(rrat,rdate,csc1,csc2,msgn,thr,wl)
%     snwdataconverter(rrat,rdate,csc1,csc2,msgn)
    close all
end