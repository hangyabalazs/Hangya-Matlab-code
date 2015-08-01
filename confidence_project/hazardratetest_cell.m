%% Dirac failure distribution

NumTrials = 10000;
ITIMin = 1;
ITIMax = 6;
ITIs = 3 * ones(1,NumTrials);   % deterministic timing

%% Failure density function

dt = 0.05;
times = ITIMin-dt:dt:ITIMax+dt;
cnts = (times(1:end-1) + times(2:end)) / 2;
ft = histc(ITIs,times);
ft = ft(1:end-1);
ft = ft / sum(ft);
figure
bar(cnts,ft)

%% Hazard rate

Ft = cumsum(ft);
Rt = 1 - Ft;
ht = (Rt(1:end-1) - Rt(2:end)) ./ (dt * Rt(1:end-1));
figure
plot(cnts(1:end-1),ht)

%% Subjective hazard rate

% Convolution of the failure density function with a time-dependent Gaussian
Phi = 0.26;
% Phi = 0.36;
intfg = zeros(size(cnts));
tau = cnts;
cntr = 1;   % counter
cnts2 = [cnts cnts(end)+dt:dt:cnts(end)+5000*dt];
for t = cnts2
    fg = ft .* exp(-(tau-t).^2/(2*Phi^2*t^2));
    intfg(cntr) = sum(fg*dt);
    cntr = cntr + 1;
end
fwt = 1 ./ (Phi * cnts2 * sqrt(2*pi)) .* intfg;
fwt = fwt / sum(fwt);

% Subjective hazard rate
Fwt = cumsum(fwt);
Rwt = 1 - Fwt;
hwt = (Rwt(1:end-1) - Rwt(2:end)) ./ (dt * Rwt(1:end-1));
figure
plot(cnts2(1:end-1),hwt)
xlim([0 7])