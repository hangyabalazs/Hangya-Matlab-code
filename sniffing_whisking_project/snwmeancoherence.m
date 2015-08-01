function snwmeancoherence
%SNWMEANCOHERENCE   Coherence between whisking and respiration.
%   SNWMEANCOHERENCE calculates coherence summary between root-mean-square
%   of whisking and respiration signal (sniffing; see SNWCOHERENCE for
%   details). Averages are calculated for individual animals over session;
%   then grand average is calculated over animals.
%
%   See also SNWMEANPHASE2 and SNWCOHERENCE.

% Import session data from Excel
global DATAPATH
tblfile = [DATAPATH 'SniffWhisk\data_selection.xlsx'];
[tbl0 tbl] = xlsread(tblfile);

% Import data
rats = unique(tbl(:,1));
numRats = length(rats);
mrCxy = [];
figure;
hold on
for iR = 1:numRats
    crat = rats(iR);
    sessinx = find(strcmp(tbl(:,1),crat));
    
    numSessions = length(sessinx);
    rCxy = [];
    for iS = 1:numSessions
        rrat = tbl{sessinx(iS),1};
        rdate = tbl{sessinx(iS),2};
        [sCxy f] = data_import(rrat,rdate);     % coherence for session
        rCxy = [rCxy; sCxy];    % all coherence functions for one rat
    end
    plot(f,nanmean(rCxy,1),'Color',rand(1,3))   % mean coherence function per rat
end



% -------------------------------------------------------------------------
function [Cxy f] = data_import(rrat,rdate)

% Define results directory
global DATAPATH
global DATADIR
resdir = [DATAPATH 'SniffWhisk\meancoherence\'];

% Load data
fnme = [DATAPATH 'SniffWhisk\coherence3\COHERENCE_' rrat '_' rdate '.mat'];
load(fnme)

% Convert to row vector
Cxy = Cxy';