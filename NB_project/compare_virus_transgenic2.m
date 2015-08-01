function compare_virus_transgenic2
%COMPARE_VIRUS_TRANSGENIC   Compare behavior of different mouse lines.
%   COMPARE_VIRUS_TRANSGENIC compares ChAT-ChR2 and ChAT-Cre mice. Two-way
%   ANOVA is calculated for go-performance, no-go-performance and go
%   reaction time with difficulty and genotype as factors.
%
%   See also PSYCH_GONOGO2.

% Load behav. data
global DATAPATH
load([DATAPATH 'NB\average_performance_newdata\PSY_RT_VARS.mat'])   % NB
GoPerformance1 = GoPerformance;
NoGoPerformance1 = NoGoPerformance;
GoReactionTime1 = GoReactionTime;
NoGoReactionTime1 = NoGoReactionTime;
load([DATAPATH 'HDB\average_performance_newdata\PSY_RT_VARS.mat'])   % HDB
GoPerformance2 = GoPerformance;
NoGoPerformance2 = NoGoPerformance;
GoReactionTime2 = GoReactionTime;
NoGoReactionTime2 = NoGoReactionTime;

% Merge
GoReactionTime = [GoReactionTime1 GoReactionTime2];
NoGoReactionTime = [NoGoReactionTime1 NoGoReactionTime2];
GoPerformance = [GoPerformance1 GoPerformance2];
NoGoPerformance = [NoGoPerformance1 NoGoPerformance2];

% Reshape
goperf = reshape(GoPerformance,3,length(GoPerformance)/3);
nogoperf = reshape(NoGoPerformance,3,length(NoGoPerformance)/3);
gort = reshape(GoReactionTime,3,length(GoReactionTime)/3);
nogort = reshape(NoGoReactionTime,3,length(NoGoReactionTime)/3);
NumMice2 = 18;
si = repmat([0; 0.5; 1],1,NumMice2);

% Mice: 13 15 18 20 23 26 27 28 29 37 38 39 40 43 45 46 71 72 77 67 70 78
% ChAT-Cre: 13 15 18 20 23 29 37 38 39 71 72 77 67 70 78
% ChAT-ChR2: 40 45 46
virusinx = [1:5 9:12 17:22];   % virus-injected
tginx = 14:16;   % transgenic
pvinx = [6:8 13];   % PV-cre
goperf = goperf(:,[virusinx tginx]);
nogoperf = nogoperf(:,[virusinx tginx]);
gort = gort(:,[virusinx tginx]);
nogort = nogort(:,[virusinx tginx]);
grp = [ones(1,15) zeros(1,3)];   % grouping variable

% ANOVA
y = reshape(goperf',NumMice2*3,1);
g1 = reshape(si',NumMice2*3,1);
g2 = [ones(15*3,1); zeros(3*3,1)];
p = anovan(y,{g1 g2})

y = reshape(nogoperf',NumMice2*3,1);
p = anovan(y,{g1 g2})

y = reshape(gort',NumMice2*3,1);
p = anovan(y,{g1 g2})

y = reshape(nogort',NumMice2*3,1);
p = anovan(y,{g1 g2})