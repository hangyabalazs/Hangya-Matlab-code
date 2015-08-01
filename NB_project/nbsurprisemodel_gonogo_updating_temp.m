function nbsurprisemodel_gonogo
%NBSURPRISEMODEL_GONOGO   HMM model go/no-go task.
%   NBSURPRISEMODEL_GONOGO implements a hidden Markov model for
%   for the auditory go/no-go inspired by the reference below.
%
%   The animal is assumed to have a statistically veridical model of the
%   task. The model is a HMM with three states: START, GO and NO-GO.
%   Leaving START is determined by the hazard function of stimulus onset,
%   with equal probabilities of GO and NO-GO states. Outputs are generated
%   from GO and NO-GO states. E.g. in the GO state a 'go' or no output is
%   generated every time step (0.01s). The animal has access only to the
%   outputs (observations) and infers state probabilities based on Bayes
%   rule. Response probability at each time stamp is defined by 
%   a*P(GO)-b*P(NO GO) where a and b are constant parameters. Reinforcement
%   surprise is defined through a hyperbolic function of the cumulative
%   number of 'go' observations before the first response (first lick
%   terminates the tone in the task).
%
%   Reference:
%   Dayan P, Yu JA (2006) Phasic norepinephrine: A neural interrupt signa;
%   for unexpected events. Network: Computation in Neural Systems
%   17(4):335-350

% we need base lick rate instead of the hazard taking care of things!!!

% Control plotting behavior
display = false;

% Parameters
detection_probability = [0.01 0.33 0.66 0.99];   % difficult to easy

% Time
dt = 0.01;   % time step
T = 0:dt:3;   % time

% Stimulus onset probability
ITId = ITIsim(200000);
ITIDistribution = hist(ITId,T);
ITIDistribution = ITIDistribution / sum(ITIDistribution) / dt;
NumTrials = 5000;
ITIs = ITIsim(NumTrials);

% Simulate trial
O = cell(1,NumTrials);   % outcomes
[tone_onset difficulty trialtype eta ...
    RT isdetected G_strength NG_strength ...
    lickgocoeff licknogocoeff] = deal(nan(1,NumTrials));
meta = mean(1-(1-detection_probability(1:4)).^(dt/0.5));   % mean detection prob. in time bin
lgdefault = 0.05;   % stady-state go lick coefficient 0.05
lngdefault = 0.03;   % stady-state no-go lick coefficient 0.03
basedefault = 0.5;
for iT = 1:NumTrials
    tone_onset(iT) = ITIs(iT);   % tone onset (foreperiod)
    difficulty(iT) = randi([1 4],1);   % 4 levels of difficulty
    trialtype(iT) = randi([1 2],1);   % trialtypes: 1 = no-go; 2 = go
    eta(iT) = 1 - (1 - detection_probability(difficulty(iT))) ^ (dt / 0.5);   % detection probability in time bin
    
    % Update lick probability (eagerness to lick)
    if iT == 1   % initial condition
        lickgocoeff(iT) = lgdefault;
        licknogocoeff(iT) = lngdefault;
    else
        switch O{iT-1}
            case 'hit'   % increase positive corr. between go prob. and lick
                lickgocoeff(iT) = lickgocoeff(iT-1) * 1.2;
                licknogocoeff(iT) = licknogocoeff(iT-1) * 0.8;
            case 'fa'    % increase negative corr. between go prob. and lick
                lickgocoeff(iT) = lickgocoeff(iT-1) * 0.8;
                licknogocoeff(iT) = licknogocoeff(iT-1) * 1.2;
            case 'cr'    % no updating
                lickgocoeff(iT) = mean([lickgocoeff(iT-1) lgdefault]);
                licknogocoeff(iT) = mean([licknogocoeff(iT-1) lngdefault]);
            case 'miss'   % no updating
                lickgocoeff(iT) = mean([lickgocoeff(iT-1) lgdefault]);
                licknogocoeff(iT) = mean([licknogocoeff(iT-1) lngdefault]);
        end
    end
    
    O{iT} = 'restart';
    while isequal(O{iT},'restart')
        [O{iT} RT(iT) isdetected(iT) G_strength(iT) NG_strength(iT)] = ...
            runtrial(T,ITIDistribution,eta(iT),tone_onset(iT),trialtype(iT),meta,...
            lickgocoeff(iT),licknogocoeff(iT),display);
    end
    disp(iT)
end

% Psychometric function & 'surprise'
K = 1;
[GoPerf NoGoPerf GoRT Go_Surprise NoGo_Surprise] = deal(nan(1,4));
for d = 1:4
    GoPerf(d) = sum(strcmp(O,'hit')&difficulty==d) / sum(trialtype==1&difficulty==d);   % response rate for go tones
    NoGoPerf(d) = sum(strcmp(O,'fa')&difficulty==d) / sum(trialtype==2&difficulty==d);   % response rate for no-go tones
    
    GoRT(d) = mean(RT(strcmp(O,'hit')&difficulty==d));
    
    Go_Surprise(d) = mean(K./(G_strength(strcmp(O,'hit')&difficulty==d)+K));   % 'surprise' based on the number of detection events
    NoGo_Surprise(d) = mean(K./(NG_strength(strcmp(O,'fa')&difficulty==d)+K));
%     Go_Surprise(d) = mean(1-G_strength(strcmp(O,'hit')&difficulty==d));   % 'surprise' based on the number of detection events
%     NoGo_Surprise(d) = mean(1-NG_strength(strcmp(O,'fa')&difficulty==d));
end

% Plot
figure
plot(GoPerf,'Color',[0 0.8 0],'Marker','o','MarkerFaceColor',[0 0.8 0],...
    'MarkerEdgeColor',[0 0.8 0],'MarkerSize',12,'LineWidth',2)
hold on
plot(NoGoPerf,'Color',[0.8 0 0],'Marker','o','MarkerFaceColor',[0.8 0 0],...
    'MarkerEdgeColor',[0.8 0 0],'MarkerSize',12,'LineWidth',2)

figure
bar(1:4,NoGo_Surprise,'FaceColor',[0.8 0 0])
hold on
bar(6:9,Go_Surprise,'FaceColor',[0 0.8 0])

figure
plot(GoRT,'k','LineWidth',2)
keyboard

% Updating
next_hit_inx = strcmp([O(2:end) NaN],'hit');
prev_hit_inx = strcmp([{NaN} O(1:end-1)],'hit');
next_fa_inx = strcmp([O(2:end) NaN],'fa');
prev_fa_inx = strcmp([{NaN} O(1:end-1)],'fa');
[GPhp NGPhp GPhm NGPhm GPfp NGPfp GPfm NGPfm] = deal(nan(1,4));
for d = 1:4
    GPhp(d) = sum(strcmp(O,'hit')&prev_hit_inx&difficulty==d) / sum(trialtype==1&prev_hit_inx&difficulty==d);
    NGPhp(d) = sum(strcmp(O,'fa')&prev_hit_inx&difficulty==d) / sum(trialtype==2&prev_hit_inx&difficulty==d);
    GPhm(d) = sum(strcmp(O,'hit')&next_hit_inx&difficulty==d) / sum(trialtype==1&next_hit_inx&difficulty==d);
    NGPhm(d) = sum(strcmp(O,'fa')&next_hit_inx&difficulty==d) / sum(trialtype==2&next_hit_inx&difficulty==d);
    GPfp(d) = sum(strcmp(O,'hit')&prev_fa_inx&difficulty==d) / sum(trialtype==1&prev_fa_inx&difficulty==d);
    NGPfp(d) = sum(strcmp(O,'fa')&prev_fa_inx&difficulty==d) / sum(trialtype==2&prev_fa_inx&difficulty==d);
    GPfm(d) = sum(strcmp(O,'hit')&next_fa_inx&difficulty==d) / sum(trialtype==1&next_fa_inx&difficulty==d);
    NGPfm(d) = sum(strcmp(O,'fa')&next_fa_inx&difficulty==d) / sum(trialtype==2&next_fa_inx&difficulty==d);
end
GoHitUpdate = GPhp - GPhm;
NoGoHitUpdate = NGPhp - NGPhm;
GoFAUpdate = GPfp - GPfm;
NoGoFAUpdate = NGPfp - NGPfm;
figure
plot(GoHitUpdate,'g')
hold on
plot(NoGoHitUpdate,'r')
figure
plot(GoFAUpdate,'g')
hold on
plot(NoGoFAUpdate,'r')
keyboard

% Hazard effect on surprise
inx_unexpected = tone_onset > 0.5 & tone_onset < 1.5;
inx_expected = ~inx_unexpected;
Go_Surprise_expected = mean(K./(G_strength(strcmp(O,'hit')&inx_expected)+K));   % 'surprise' based on the number of detection events
Go_Surprise_unexpected = mean(K./(G_strength(strcmp(O,'hit')&inx_unexpected)+K));   % 'surprise' based on the number of detection events
Go_Surprise_expected = mean(1-G_strength(strcmp(O,'hit')&inx_expected));   % 'surprise' based on GO probability
Go_Surprise_unexpected = mean(1-G_strength(strcmp(O,'hit')&inx_unexpected));   % 'surprise' based on GO probability


% -------------------------------------------------------------------------
function [O RT isdetected G_strength NG_strength] = runtrial(T,ITIDistribution,eta,tone_onset,trialtype,meta,lickgocoeff,licknogocoeff,display)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Hidden Markov model                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k = 0;
dt = T(2) - T(1);   % time step
[S X] = deal(nan(size(T)));   % STATES, OUTPUTS
for t = T
    k = k + 1;
    
    if t < tone_onset
        S(k) = 0;   % STATE: 0 = start
    elseif t >= tone_onset && t <= tone_onset + 0.5   % tone durstion: 0.5s
        S(k) = trialtype;   % STATE: 1 = no-go; 2 = go
    else
        S(k) = 3;   % STATE: 3 = terminate
    end
    
    switch S(k)
        case 0   % start STATE
            X(k) = 0;   % OUTPUT: 0 = no output
        case 1   % go STATE
            r = rand;
            if r < eta
                X(k) = 1;   % OUTPUT: 1 = G
            else
                X(k) = 0;   % OUTPUT: 0 = no output
            end
        case 2  % no-go STATE
            r = rand;
            if r < eta
                X(k) = 2;   % OUTPUT: 2 = NG
            else
                X(k) = 0;   % OUTPUT: 0 = no output
            end
        case 3  % terminate STATE
            X(k) = 0;
    end
end
isdetected = any(X>0);

% Plot
I1 = cumsum(X==1);
I2 = cumsum(X==2);
if display
    figure
    plot(T,S,'k')   % STATE
    
    figure
    plot(T,X,'ko')   % OUTPUT
    
    figure
    plot(T,I1.*dt,'g');   % cumulative observation of G
    hold on
    plot(T,I2.*dt,'r');    % cumulative observation of NG
    legend({'NG' 'G'})
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Exact probabilistic inference                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p_Xt_St = [1 0 0; 1-meta meta 0; 1-meta 0 meta];   % p(x_t|s_t) output probability matrix
% row 1: s_t = 0
% row 2: s_t = 1
% row 3: s_t = 2
% column 1: x_t = 0
% column 2: x_t = 1
% column 3: x_t = 2

% Inference
k = 1;
[Ps Pg Png q] = deal(nan(size(T)));   % state probabilities based on output observations
Ps(1) = 1;   % prob. of start STATE: 1 in the first step
Pg(1) = 0;   % prob. of go STATE: 0 in the first step
Png(1) = 0;   % prob. of no-go STATE: 0 in the first step
q(1) = 0;   % p(tone(k)|no tone(1,...,k-1))
for t = T(2:end)
    k = k + 1;
    q(k) = ITIDistribution(k) / sum(ITIDistribution(k:end));   % hazard
    Ps(k) = p_Xt_St(1,X(k)+1) * (1-q(k)) * Ps(k-1);
    Pg(k) = p_Xt_St(2,X(k)+1) * (0.5 * q(k) * Ps(k-1) + Pg(k-1));
    Png(k) = p_Xt_St(3,X(k)+1) * (0.5 * q(k) * Ps(k-1) + Png(k-1));
end

% Normalize probabilities
C = Ps + Pg + Png;
Ps = Ps ./ C;
Pg = Pg ./ C;
Png = Png ./ C;

% Plot
if display
    figure
    subplot(3,1,1)
    plot(T,Ps,'k');    % P(start)
    ylim([0 1])
    legend({'\itP(start)'})
    subplot(3,1,2)
    plot(T,Pg,'k');    % P(go)
    ylim([0 1])
    legend({'\itP(go)'})
    subplot(3,1,3)
    plot(T,Png,'k');    % P(no-go)
    ylim([0 1])
    legend({'\itP(no-go)'})
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                  Lick                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Animal response
base_lick_rate = 1 - (1 - 0.2) ^ (dt / 0.5);
P_lick = 0 * base_lick_rate + lickgocoeff * Pg - licknogocoeff * Png;   % lick probability (1; 0.035; 0.032)
% P_lick = 1 * base_lick_rate + lickgocoeff * Pg - licknogocoeff * Png;
P_lick(P_lick<0) = 0;
r = rand(size(T));
R = r < P_lick;   % lick
alllickinx = find(R);

% Reaction time
lick_TS = T(R);
lickinx = find(lick_TS>tone_onset&lick_TS<tone_onset+0.5,1,'first');
first_lick_TS = lick_TS(lickinx);
RT = first_lick_TS - tone_onset;
if isempty(RT)
    RT = NaN;
end

% Stimulus strength
if ~isnan(RT)
    G_strength = I1(alllickinx(lickinx));   % G 'strength': cumulative number of 'go' outputs
    NG_strength = I2(alllickinx(lickinx));   % NG 'strength': cumulative number of 'no go' outputs
%     G_strength = Pg(alllickinx(lickinx));
%     NG_strength = Png(alllickinx(lickinx));
else   % no lick: miss or correct rejection
    G_strength = I1(end);   % G 'strength': cumulative number of 'go' outputs
    NG_strength = I2(end);   % NG 'strength': cumulative number of 'no go' outputs
end

% Plot
if display
    figure
    subplot(2,1,1)
    plot(T,P_lick,'k')
    subplot(2,1,2)
%     plot(T,R,'k')
    NumLicks = length(lick_TS);
    line([lick_TS; lick_TS],[zeros(1,NumLicks); ones(1,NumLicks)],'Color','k')
    xlim([0 3])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Trial outcome                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Trial outcome
if any(R&S==0)
    O = 'restart';   % restart trial
else
    if trialtype == 1   % go STATE
        if any(R&S==1)   % animal response
            O = 'hit';   % hit
        else   % no response
            O = 'miss';   % miss
        end
    else   % no-go STATE
        if any(R&S==2)   % animal response
            O = 'fa';   % false alarm
        else   % no response
            O = 'cr';   %  correct rejection
        end
    end
end
if display
    disp(O)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               ACh signal                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ACh = Pg / 0.2;   % P(target|X) / P(target)
% 
% % Plot
% if display
%     figure
%     plot(T,ACh,'k')
%     ylim([0 5])
% end

% -------------------------------------------------------------------------
function ITIs = ITIsim(NumTrials)

% Parameters of the foreperiod distribution
ITIMin = 0.1;   % parameters of the uniform distribution
ITIMax = 3;
mng1 = 0.3;   % parameters for the Gaussians
mng2 = 2;
sdg = 0.15;
pmx1 = 0.35;   % mixing probabilities
pmx2 = 0.35;
pmx3 = 1 - pmx1 - pmx2;

% Simulate
% Sample the failure distribution
ITIs1 = random('Normal',mng1,sdg,1,NumTrials);
while any(ITIs1>ITIMax) || any(ITIs1<ITIMin)
    inx = ITIs1 > ITIMax  | ITIs1 < ITIMin;
    ITIs1(inx) = random('Normal',mng1,sdg,1,sum(inx));
end
ITIs2 = random('Normal',mng2,sdg,1,NumTrials);
while any(ITIs2>ITIMax) || any(ITIs2<ITIMin)
    inx = ITIs2 > ITIMax  | ITIs2 < ITIMin;
    ITIs2(inx) = random('Normal',mng2,sdg,1,sum(inx));
end
ITIs3 = random('Uniform',ITIMin,ITIMax,1,NumTrials);
prr = rand(1,NumTrials);
rr = zeros(3,NumTrials);
rr(1,prr<pmx1) = 1;
rr(2,prr>=pmx1&prr<(pmx1+pmx2)) = 1;
rr(3,prr>=(pmx1+pmx2)) = 1;
ITIs = rr(1,:) .* ITIs1 + rr(2,:) .* ITIs2 + rr(3,:) .* ITIs3;