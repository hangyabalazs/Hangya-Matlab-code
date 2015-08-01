function snwmeancoherence2
%SNWMEANCOHERENCE2   Coherence between whisking and respiration.
%   SNWMEANCOHERENCE2 calculates coherence summary between root-mean-square
%   of whisking and respiration signal (sniffing; see SNWCOHERENCE for
%   details). Averages are calculated for individual animals over session;
%   then grand average is calculated over animals.
%
%   See also SNWMEANPHASE2 and SNWCOHERENCE.

% Import session data from Excel
global DATAPATH
tblfile = [DATAPATH 'SniffWhisk\data_selection_coherence.xlsx'];
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
        e1 = tbl0(sessinx(iS),4);
        e2 = tbl0(sessinx(iS),5);
        [sCxy f] = data_import(rrat,rdate,e1,e2);     % coherence for session
        rCxy = [rCxy; sCxy];    % all coherence functions for one rat
    end
    plot(f,nanmean(rCxy,1),'Color',rand(1,3))   % mean coherence function per rat
    if f(end) == 25
        mrCxy = [mrCxy; nanmean(rCxy(:,1:2:end),1)];
        F = f(1:2:end);
    elseif f(end) == 50
        mrCxy = [mrCxy; nanmean(rCxy(:,1:65),1)];
    end
end

keyboard

% Plot
figure
plot(F,nanmean(mrCxy))



% -------------------------------------------------------------------------
function [Cxy f] = data_import(rrat,rdate,e1,e2)

% Define results directory
global DATAPATH
global DATADIR
resdir = [DATAPATH 'SniffWhisk\meancoherence\'];

% Load data
if isnan(e1)
    e1 = 'start';
end
if isnan(e2)
    e2 = 'end';
end
fnme = [DATAPATH 'SniffWhisk\coherence_selected4\COHERENCE_' rrat '_' ...
    rdate '_' num2str(e1) '_' num2str(e2) '.mat'];
load(fnme)

% Convert to row vector
Cxy = Cxy';