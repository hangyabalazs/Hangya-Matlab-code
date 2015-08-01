function nbsurprisemodel_gonogo
%NBSURPRISEMODEL_DAYAN   HDD model of NE activity.
%   NBSURPRISEMODEL_DAYAN implements the hidden Markov model for
%   norepinephrine release presented in the reference below.
%
%   Reference:
%   Daya P, Yu JA (2006) Phasic norepinephrine: A neural interrupt signa;
%   for unexpected events. Network: Computation in Neural Systems
%   17(4):335-350

% Parameters
E = 0.1;   % expectancy parameter
Eta = [0.7 0.8 0.9 1];   % difficulty (difficult to easy)

% Time
dt = 0.01;   % time step
T = 0:dt:3;   % time

% Stimulus onset probability
ITId = ITIsim(20000);
ITIDistribution = hist(ITId,T);
ITIDistribution = ITIDsitribution / sum(ITIDistribution) / dt;
NumTrials = 600;
ITIs = ITIsim(NumTrials);

% Simulate trial
for iT = 1:NumTrials
    tone_onset = ITIs(iT);   % tone onset (foreperiod)
    difficulty = randi([1 4],1);   % 4 levels of difficulty
    trialtype = randi([1 2],1);   % trialtypes: 1 = no-go; 2 = go
    runtrial(T,ITIDistribution,Eta,tone_onset,difficulty,trialtype)
end

% -------------------------------------------------------------------------
function runtrial(T,q,Eta,tone_onset,difficulty)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Hidden Markov model                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k = 0;
[S X] = deal(nan(size(T)));   % STATES, OUTPUTS
for t = T
    k = k + 1;
    
    if t < tone_onet
        S(k) = 0;   % STATE: 0 = start
    elseif t >= tone_onset && t <= tone_onset + 0.5   % tone durstion: 0.5s
        if k == 1 || S(k-1) == 0
            S(k) = trialtype;   % STATE: 1 = no-go; 2 = go
        end
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
    end
end

% Plot
figure
plot(T,S,'k')   % STATE

figure
plot(T,X,'ko')   % OUTPUT

figure
plot(T,cumsum(X==1).*dt,'g');   % cumulative observation of G
hold on
plot(T,cumsum(X==2).*dt,'r');    % cumulative observation of NG
ylim([0 T(end)-10*dt])
legend({'NG' 'G'})

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Exact probabilistic inference                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p_Xt_St = [1 0 0; 1-eta eta 0; 1-eta 0 eta];   % p(x_t|s_t) output probability matrix
% row 1: s_t = 0
% row 2: s_t = 1
% row 3: s_t = 2
% column 1: x_t = 0
% column 2: x_t = 1
% column 3: x_t = 2

% Inference
k = 1;
[Ps Pg Png] = deal(nan(size(T)));   % state probabilities based on output observations
Ps(1) = 1;   % prob. of start STATE: 1 in the first step
Pg(1) = 0;   % prob. of go STATE: 0 in the first step
Png(1) = 0;   % prob. of no-go STATE: 0 in the first step
for t = T(2:end)
    k = k + 1;
    Ps(k) = p_Xt_St(1,X(k)+1) * (Ps(k-1) * P_notone + (1-q(k)) * Ps(k-1) * P_tone);
    Pg(k) = p_Xt_St(2,X(k)+1) * (0.5 * q(k) * Ps(k-1) + Pg(k-1));
    Png(k) = p_Xt_St(3,X(k)+1) * (0.5 * q(k) * Ps(k-1) + Png(k-1));
end

% Normalize probabilities
C = Ps + Pg + Png;
Ps = Ps ./ C;
Pg = Pg ./ C;
Png = Png ./ C;

% Plot
figure
subplot(3,1,1)
plot(T,Ps,'k');    % P(start)
ylim([0 1])
legend({'\itP(start)'})
subplot(3,1,2)
plot(T,Pd,'k');    % P(distract)
ylim([0 1])
legend({'\itP(distract)'})
subplot(3,1,3)
plot(T,Pt,'k');    % P(target)
ylim([0 1])
legend({'\itP(target)'})

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Trial outcome                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Target probability reaches threshold
if any(Pt>0.95)
    response_time1 = T(find(Pt>0.95,1,'first'));   % upper limit: response
else
    response_time1 = Inf;
end
if any(Pt<0.01&Ps==0)
    response_time0 = T(find(Pt<0.01&Ps==0,1,'first'));   % lower limit: no response
else
    response_time0 = Inf;
end

% Animal response
if response_time1 < response_time0
    R = 1;   % animal response
else
    R = 0;   % no response
end

% Trial outcome
if S(end) == 1   % distract STATE 
    if R == 1   % animal response
        O = 'fa';   % false alarm
    else   % no response
        O = 'cr';   % correct rejection
    end
else   % target STATE
    if R == 1   % animal response
        O = 'hit';   % hit
    else   % no response
        O = 'miss';   %  miss
    end
end
disp(O)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               ACh signal                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NE = Pt / 0.2;   % P(target|X) / P(target)

% Plot
figure
plot(T,NE,'k')
ylim([0 5])

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

% Time
if nargin < 2
    dt = 0.01;     % resolution
end
t = ITIMin:dt:ITIMax;   % time

% Simulate
if nargout == 4
    
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
end