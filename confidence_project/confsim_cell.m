%% percept distribution

% Generate series (distribution) of percepts
N = 1000000;
Rd = randi([0 1],1,N);
L = Rd + 1;  % left click rate; 1 or 2
R = 2 - Rd;  % right click rate; 2 or 1
L_hat = 1 * randn(1,N) + L;  % left click rate
R_hat = 1 * randn(1,N) + R ;  % right click rate
D_hat = L_hat - R_hat;  % perceived evidence ('percept')
Perf = sum((D_hat>0&L>R)|(D_hat<0&L<R)) / N;   % average performance

%% confidence

% Discretize percept
mn = min(D_hat);
mx = max(D_hat);
bno = 100;  % bin number: d_hat has to be discretized
D_hat_disc = round((D_hat-mn)/(mx-mn)*bno);

% Choice: sign of percept
% Choice = double(D_hat>0);   % choice distribution: 1=L, 0=R

%% simulate trials

% Trials
NumTrials = 10000;  % number of trials
[conf outcome] = deal(nan(1,NumTrials));
for iT = 1:NumTrials   % loop through trials
    
    % Draw click rates; uniform prior
    rd = randi([0 1],1,1);
    l = rd + 1;  % left click rate; 1 or 2
    r = 2 - rd;  % right click rate; 2 or 1
    
    % Perceived evidence (percept) per trial
    l_hat = 1 * randn(1,1) + l;  % left click rate
    r_hat = 1 * randn(1,1) + r ;  % right click rate
    d_hat = l_hat - r_hat;  % perceived evidence ('percept')
    
    % Choice: sign of percept
    choice = double(d_hat>0);   % choice: 1=L, 0=R
    
    % Confidence per trial
    d_hat_inx = round((d_hat-mn)/(mx-mn)*bno);  % bin index for discretized d_hat
    if choice   % left side choice
        P_H_L_d_hat = sum(L>R&D_hat_disc==d_hat_inx) / N;   % P(H_L,d_hat)
        P_d_hat = sum(D_hat_disc==d_hat_inx) / N;   % P(d_hat)
        conf(iT) = P_H_L_d_hat / P_d_hat;   % confidence per trial
    else   % right side choice
        P_H_R_d_hat = sum(R>L&D_hat_disc==d_hat_inx) / N;   % P(H_R,d_hat)
        P_d_hat = sum(D_hat_disc==d_hat_inx) / N;   % P(d_hat)
        conf(iT) = P_H_R_d_hat / P_d_hat;   % confidence per trial
    end
    
    % Outcome per trial
    if (choice & l>r) | (~choice & r>l) %#ok<OR2,AND2>
        outcome(iT) = 1;  % correct trial
    else
        outcome(iT) = 0;   % incorrect trial
    end
    
end

%% accuracy (percent correct) vs confidence)

% Accuracy for a given confidence level: A(c)
A = [];
mnc = 0.5;
mxc = 1;
bnoc = 5;  % bin number: discretize confidence
conf_disc = round((conf-mnc)/(mxc-mnc)*bnoc);   % discrete confidence values
conf_levels = unique(conf_disc);  % confidence levels
conf_levels = conf_levels(conf_levels>0);  % in a few cases conf. < 0.5
nconf = bnoc;   % number of different confidence levels
[A_c mconf] = deal(nan(1,nconf));
for iC = 1:nconf   % calculate accuracy for each conf. level
    P_c_correct = sum(conf_disc==conf_levels(iC)&outcome==1) / NumTrials;  % P(c, correct)
    P_c = sum(conf_disc==conf_levels(iC)) / NumTrials;   % P(c)
    A_c(iC) = P_c_correct / P_c;  % accuracy: A(c)
    
    mconf(iC) = mean(conf(conf_disc==conf_levels(iC)));  % confidence value corresponding to each level
end

% Plot
figure
plot(conf_levels,A_c)  % by level

figure
plot(mconf,A_c)  % by value