function rastersim
%RASTERSIM   Adaptive PSTH calculation.
%
%   See also SSKERNEL.

% Probabilities
% p = [5 5 80 2 3 2 6 6 7 6 8 10 10 9 12 11 12 13 14 15 14 13 12 11 10 10 ...
%     9 8 7 7 6 6 5 5 5 5 5 5 5 5 5 5 5 5 5] / 100;
p = [2 3 2 6 6 7 6 8 10 10 9 12 11 12 13 14 15 14 13 12 11 10 10 ...
    9 8 7 7 6 6 5 5 5 5 5 5 5 5 5 5 5 5 5] / 100;
p2 = repmat(p,30,1) / 10;
p = [ones(1,10)*0.05 0.75 p2(:)'];
tl = length(p);

% Simulate trials
tno = 30;   % number of trials
spt = zeros(tno,tl);
for k = 1:tno
    spt(k,:) = double(rand(1,tl)<p);
end
figure
imagesc(spt)
colormap('bone')

% PSTH
psth = sum(spt);
figure
plot(psth)

% Convolution with Gaussian kernel
% wbh = gausswin(20);
% gvd = zeros(tno,tl);
% for k = 1:tno   % convolve trial-wise
%     pgvd = conv(spt(k,:),wbh);
%     lag = (length(pgvd)-tl+1) / 2;
%     gvd(k,:) = pgvd(lag:end-lag);
% end
% figure
% imagesc(gvd)
% 
% psth_gconv = sum(gvd);
% figure
% plot(psth_gconv)

% Convolution with variable Gaussian kernel
agvd = zeros(tno,tl);
prob = sum(spt) / tno;
figure
hold on
for k = 1:tno   % convolve trial-wise
    spks = find(spt(k,:));
    spno = length(spks);
    for t = 1:spno
        spi = spks(t);
        tspt = zeros(1,tl);
        tspt(spi) = 1;
        wbh = gausswin(20,2+prob(spi)*10);
        wbh = wbh / sum(wbh);
        plot(wbh,'Color',rand(1,3))
        ptc = conv(tspt,wbh);
        lag = (length(ptc)-tl+1) / 2;
        ptc = ptc(lag:end-lag);
        agvd(k,:) = agvd(k,:) + ptc;
    end
end
figure
imagesc(agvd)

psth_aconv = sum(agvd);
figure
plot(psth_aconv)

% Convolution with variable exponential kernel
aevd = nan(1,tl);
alpha = nan(1,tl);
spts = sum(spt);    % convolution of the sum should be the same as convolving trialwise and sum up (might be worth to check)
alpha(1) = 0;
aevd(1) = spts(1);
for t = 2:tl   % convolve time-wise
    inx = t-min(9,t-1):t;
    prbs = sum(spts(inx));
    prob = prbs / (length(inx) * tno);
    if prob > 1
        keyboard
    end
    alpha(t) = 1-(1-prob)^6;
    aevd(t) = (1 - alpha(t)) * aevd(t-1) + alpha(t) * spts(t);
end
figure
plot(aevd)

% Shimazaki & Shinomoto
allspks = [];
for k = 1:tno
    spks = find(spt(k,:));
    allspks = [allspks spks];
end
ts = sort(allspks);     % merged spike train
[optW, C, W] = sskernel(ts);

% Shimazaki & Shinomoto
% Step1
allspks =[];
for k = 1:tno
    spks = find(spt(k,:));
    allspks = [allspks spks];
end
ts = sort(allspks);     % merged spike train
aspno = length(ts);     % number of all spikes

% Step2
W = 20;

Gamma = 0.15:0.02:0.85;
lG = length(Gamma);
C_n = zeros(1,lG);
next_gamma = 1;
for gamma = Gamma
    disp(gamma)
    W_gamma = zeros(1,tl);
    w_bar_gamma = zeros(1,tl);
    for t = 1:tl
        disp(t)
        T =0:0.01:0.4;
        lT = length(T);
        C_W_t = zeros(1,lT);
        next_w = 1;
        for w = T
            A2 = 2 * w^2 * (w^2 + 2 * W^2);
            psi_w_W_t = zeros(aspno,aspno);
            pB = zeros(aspno,aspno);
            for i = 1:aspno
                for j = 1:aspno
                    ti = ts(i);
                    tj = ts(j);
                    A1 = ((t - ti)^2 + (t - tj)^2) * w^2 + (ti - tj)^2 * W^2;
                    A = A1 / A2;
                    psi_w_W_t(i,j) = (2 * pi * w * sqrt(w^2+2*W^2))^(-1) * exp(-A);

                    pB(i,j) = funk(ti-tj,w) * funrho(W,ti);
                end
            end
            pB2 = pB * abs(1-eye(aspno));
            B = sum(pB2(:));
            C_W_t(next_w) = (1 / tno^2) * sum(psi_w_W_t(:)) - (2 / tno^2) * B;
            next_w = next_w + 1;
        end
        w_star = T(C_W_t==nanmin(C_W_t));

        %     gamma = 0.8;
        W_gamma(t) = w_star / gamma;
        w_bar_gamma(t) = w_star;
    end

    pw_gamma = zeros(tl,tl);
    w_gamma = zeros(1,tl);
    for t = 1:tl
        for s = 1:tl
            pw_gamma(s,t) = funrho(W_gamma(s),s) .* w_bar_gamma(s);
        end
        w_gamma(t) = sum(pw_gamma,1);
    end

    lambda = zeros(1,tl);
    for t = 1:tl
        pl = zeros(1,aspno);
        for i = 1:aspno
            pl(i) = funk(t-ti,w_gamma(t));
        end
        lambda(t) = sum(pl);
    end
    C = sum(lambda.^2);
    pD = zeros(aspno,aspno);
    for i = 1:aspno
        for j = 1:aspno
            ti = ts(i);
            tj = ts(j);
            pD(i,j) = funk(ti-tj,w_gamma(tj));
        end
    end
    pD2 = pD * abs(1-eye(aspno));
    D = 2 / tno^2 * sum(pD2(:));
    C_n(next_gamma) = C - D;
    next_gamma = next_gamma + 1;
end

gamma_star = Gamma(C_n==min(C_n));

% Step 4b: repeat step 2 with gamma_star
W_gamma = zeros(1,tl);
w_bar_gamma = zeros(1,tl);
for t = 1:tl
    T =0:0.01:0.4;
    lT = length(T);
    C_W_t = zeros(1,lT);
    next_w = 1;
    for w = T
        A2 = 2 * w^2 * (w^2 + 2 * W^2);
        psi_w_W_t = zeros(aspno,aspno);
        pB = zeros(aspno,aspno);
        for i = 1:aspno
            for j = 1:aspno
                ti = ts(i);
                tj = ts(j);
                A1 = ((t - ti)^2 + (t - tj)^2) * w^2 + (ti - tj)^2 * W^2;
                A = A1 / A2;
                psi_w_W_t(i,j) = (2 * pi * w * sqrt(w^2+2*W^2))^(-1) * exp(-A);

                pB(i,j) = funk(ti-tj,w) * funrho(W,ti);
            end
        end
        pB2 = pB * abs(1-eye(aspno));
        B = sum(pB2(:));
        C_W_t(next_w) = (1 / tno^2) * sum(psi_w_W_t(:)) - (2 / n^2) * B;
        next_w = next_w + 1;
    end
    w_star = w(C_W_t==min(C_W_t));

    gamma = gamma_star;
    W_gamma(t) = w_star / gamma;
    w_bar_gamma(t) = w_star;
end

pw_gamma = zeros(tl,tl);
w_gamma_star = zeros(1,tl);
for t = 1:tl
    for s = 1:tl
        pw_gamma(s,t) = funrho(W_gamma(s),s) .* w_bar_gamma(s);
    end
    w_gamma_star(t) = sum(pw_gamma,1);      % optimal bandwidth
end

% Convolution using w_gamma_star
agvd = zeros(tno,tl);
prob = sum(spt) / tno;
for k = 1:tno   % convolve trial-wise
    spks = find(spt(k,:));
    spno = length(spks);
    for t = 1:spno
        spi = spks(t);
        tspt = zeros(1,tl);
        tspt(spi) = 1;
        wbh = zeros(1,tl);
        for s = 1:tl
            wbh(s) = funk(s,w_gamma_star(t));
        end
        wbh = wbh / sum(wbh);
        ptc = conv(tspt,wbh);
        lag = (length(ptc)-tl+1) / 2;
        ptc = ptc(lag:end-lag);
        agvd(k,:) = agvd(k,:) + ptc;
    end
end
figure
imagesc(agvd)

psth_agconv = sum(agvd);
figure
plot(psth_agconv)
keyboard


% -------------------------------------------------------------------------
function r = funk(s,w)

r = 1 / (sqrt(2*pi) * w) * exp(-s^2/(2*w^2));

% -------------------------------------------------------------------------
function r = funrho(W,ti)

if ti < W / 2
    pwn = gausswin(W);
    wn = pwn / sum(pwn);
    r = wn(W/2+ti);
else
    r = 0;
end