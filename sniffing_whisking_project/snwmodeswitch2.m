function snwmodeswitch2
%SNWMODESWITCH2   Coupling mode statistics.
%   SNWMODESWITCH2 compares sniffing frequency and animal speed during
%   differnet modes of sniffing-whisking phase locking (see SNWSNIFFPLOT2).
%   It calculates statistics restricted to mode switches as well.
%
%   See also SNWSNIFFPLOT2 and SNWMODESWITCH.

% Sampling rate
sr = 1000;

% Load all segments
global DATADIR
fn = [DATADIR 'HSW\allsegments_speed_filtered.mat'];
load(fn)
fn = [DATADIR 'HSW\switches.mat'];
load(fn)

% Mode indices
gr = [allsegments.group];   % group codes corresponding to all segments
gr05 = find(gr==0.5);   % indices for 1:2 mode
gr1 = find(gr==1);   % indices for 1:1 mode
gr2 = find(gr==2);   % indices for 2:1 mode

% BAD!!!! - switch can be only whithin a session
% grc = gr(1:end-1);
% grnext = gr(2:end);
% sw1to2 = find(grc==1&grnext==2);   % mode switch from 1:1 to 2:1

% Mode frequency statistics
fr = [allsegments.frequency];   % frequencies corresponding to all segments
figure
boxplot([fr(gr05) fr(gr1) fr(gr2)],[zeros(size(gr05)) ones(size(gr1)) ones(size(gr2))*2],'labels',[{'1:2'} {'1:1'} {'2:1'}]);

% Comparison of sniffing frequency before and after the switches
boxstat(fr(sw1to2),fr(sw1to2+1),'before switch','after switch')

% Mode speed statistics
sp = cellfun(@mean,{allsegments.speed});   % frequencies corresponding to all segments
figure
boxplot([sp(gr05) sp(gr1) sp(gr2)],[zeros(size(gr05)) ones(size(gr1)) ones(size(gr2))*2],'labels',[{'1:2'} {'1:1'} {'2:1'}]);

% Statistics
sp05 = sp(gr05);
sp05(isnan(sp05)) = [];
sp1 = sp(gr1);
sp1(isnan(sp1)) = [];
sp2 = sp(gr2);
sp2(isnan(sp2)) = [];
ranksum(sp05,sp1)
ranksum(sp1,sp2)
ranksum(sp05,sp2)