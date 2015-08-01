function compare_virus_transgenic
%COMPARE_VIRUS_TRANSGENIC   Compare behavior of different mouse lines.
%   COMPARE_VIRUS_TRANSGENIC compares ChAT-ChR2 and ChAT-Cre mice. Two-way
%   ANOVA is calculated for go-performance, no-go-performance and go
%   reaction time with difficulty and genotype as factors.
%
%   See also PSYCH_GONOGO2.

% Load behav. data
load('C:\Balazs\_analysis\NB\average_performance\virus_vs_transgenic\vars.mat')

% Reshape variables
NumMice = length(amGoPerformance) / 3;
goperf = reshape(amGoPerformance,3,NumMice)';
nogoperf = reshape(amNoGoPerformance,3,NumMice)';
gort = reshape(amGoReactionTime,3,NumMice)';
nogort = reshape(amNoGoReactionTime,3,NumMice)';
si = reshape(Ib,3,NumMice)';

% Mice
mice2 = mice([2 4:6 8:length(mice)]);   % mice included
virusinx = [1:5 9:12];   % virus-injected
tginx = 14:16;   % transgenic
pvinx = [6:8 13];   % PV-cre
goperf(pvinx,:) = [];   % drop PV-cre mice
nogoperf(pvinx,:) = [];
gort(pvinx,:) = [];
nogort(pvinx,:) = [];
si(pvinx,:) = [];
mice2(pvinx) = [];
grp = [ones(1,9) zeros(1,3)];   % grouping variable

% ANOVA
NumMice2 = length(mice2);
y = reshape(goperf',NumMice2*3,1);
g1 = reshape(si',NumMice2*3,1);
g2 = [ones(9*3,1); zeros(3*3,1)];
p = anovan(y,{g1 g2})

y = reshape(nogoperf',NumMice2*3,1);
p = anovan(y,{g1 g2})

y = reshape(gort',NumMice2*3,1);
p = anovan(y,{g1 g2})

y = reshape(nogort',NumMice2*3,1);
p = anovan(y,{g1 g2})