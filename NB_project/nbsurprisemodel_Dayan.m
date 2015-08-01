%NBSURPRISEMODEL_DAYAN   HMM model of NE activity.
%   NBSURPRISEMODEL_DAYAN implements the hidden Markov model for
%   norepinephrine release presented in the reference below.
%
%   Reference:
%   Dayan P, Yu JA (2006) Phasic norepinephrine: A neural interrupt signa;
%   for unexpected events. Network: Computation in Neural Systems
%   17(4):335-350

%% parameters

dt = 0.01;   % time step
T = 0:dt:0.3;   % time
q = zeros(size(T));
q(6:10) = 1./(11-(6:10));   % temporal hazard fcn. (state transition prob.)
eta = 0.65;   % difficulty

%% simulate trial

% Simulate
k = 0;
[S X] = deal(nan(size(T)));   % STATES, OUTPUTS
for t = T
    k = k + 1;
    
    if k == 1 || S(k-1) == 0
        r = rand;
        if r > 1 - q(k)
            r = rand;
            if r < 0.8
                S(k) = 1;   % STATE: 1 = distract
            else
                S(k) = 2;   % STATE: 2 = target
            end
        else
            S(k) = 0;   % STATE: 0 = start
        end
    else
        S(k) = S(k-1);
    end
    
    switch S(k)
        case 0   % start STATE
            X(k) = 0;   % no OUTPUT
        case 1   % distract STATE
            r = rand;
            if r < eta
                X(k) = 1;   % OUTPUT: 1 = D
            else
                X(k) = 2;   % OUTPUT: 2 = T
            end
        case 2  % target STATE
            r = rand;
            if r < eta
                X(k) = 2;   % OUTPUT: 2 = T
            else
                X(k) = 1;   % OUTPUT: 1 = D
            end
    end
end

% Plot
figure
plot(T,S,'k')   % STATE

figure
plot(T,X,'ko')   % OUTPUT

figure
plot(T,cumsum(X==2).*dt,'k');    % cumulative observation of T
hold on
plot(T,cumsum(X==1).*dt,'k:');   % cumulative observation of D
ylim([0 T(end)-10*dt])
legend({'T' 'D'})

%% exact probabilistic inference

p_Xt_St = [1 0 0; 0 eta 1-eta; 0 1-eta eta];   % p(x_t|s_t) output probability matrix
% row 1: s_t = 0
% row 2: s_t = 1
% row 3: s_t = 2
% column 1: x_t = 0
% column 2: x_t = 1
% column 3: x_t = 2

% Inference
k = 1;
[Ps Pd Pt] = deal(nan(size(T)));   % state probabilities based on output observations
Ps(1) = 1;   % prob. of start STATE: 1 in the first step
Pd(1) = 0;   % prob. of distract STATE: 0 in the first step
Pt(1) = 0;   % prob. of target STATE: 0 in the first step
for t = T(2:end)
    k = k + 1;
    Ps(k) = p_Xt_St(1,X(k)+1) * (1-q(k)) * Ps(k-1);
    Pd(k) = p_Xt_St(2,X(k)+1) * (0.8 * q(k) * Ps(k-1) + Pd(k-1));
    Pt(k) = p_Xt_St(3,X(k)+1) * (0.2 * q(k) * Ps(k-1) + Pt(k-1));
end

% Normalize probabilities
C = Ps + Pd + Pt;
Ps = Ps ./ C;
Pd = Pd ./ C;
Pt = Pt ./ C;

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

%% trial outcome

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

%% NE signal

NE = Pt / 0.2;   % P(target|X) / P(target)

% Plot
figure
plot(T,NE,'k')
ylim([0 5])