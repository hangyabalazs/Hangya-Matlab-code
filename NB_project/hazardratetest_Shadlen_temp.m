%% Shadlen bimodal failure distribution

NumTrials = 10000;
ITIs1 = random('Rayleigh',1/6,1,NumTrials) + 0.1;
ITIs2 = random('Rayleigh',1/sqrt(30),1,NumTrials) + 1.75;
rr = round(rand(1,NumTrials));
ITIs = rr .* ITIs1 + (1 - rr) .* ITIs2;
ITIMin = min(ITIs);
ITIMax = max(ITIs);

figure
hist(ITIs,50)

%% My bimodal failure distribution

NumTrials = 10000;
ITIMin = 0.1;
ITIMax = 50;
mng1 = 0.3;   % parameters for the Gaussians
mng2 = 2;
sdg = 0.15;
pmx1 = 0.35;   % mixing probabilities
pmx2 = 0.35;
pmx3 = 1 - pmx1 - pmx2;

ITIs1 = random('Normal',mng1,sdg,1,NumTrials);
while any(ITIs1>ITIMax) | any(ITIs1<ITIMin)
    inx = ITIs1 > ITIMax  | ITIs1 < ITIMin;
    ITIs1(inx) = random('Normal',mng1,sdg,1,sum(inx));
end
ITIs2 = random('Normal',mng2,sdg,1,NumTrials);
while any(ITIs2>ITIMax) | any(ITIs2<ITIMin)
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
ITIMin = min(ITIs);
ITIMax = max(ITIs);

figure
hist(ITIs,50)

%% My unimodal failure distribution

NumTrials = 10000;
ITIMin = 0.1;
ITIMax = 3;
mng = 1.4;   % parameters for the Gaussian
sdg = 0.25;

ITIs = random('Normal',mng,sdg,1,NumTrials);
while any(ITIs>ITIMax) | any(ITIs<ITIMin)
    inx = ITIs > ITIMax  | ITIs < ITIMin;
    ITIs(inx) = random('Normal',mng,sdg,1,sum(inx));
end

figure
hist(ITIs,50)

%% My unimodal failure distribution #2

NumTrials = 10000;
ITIMin = 0.1;
ITIMax = 3;
mng = 1.4;   % parameters for the Gaussian
sdg = 0.25;
pmx1 = 0.65;   % mixing probabilities
pmx2 = 1 - pmx1;

ITIs1 = random('Normal',mng,sdg,1,NumTrials);
while any(ITIs1>ITIMax) | any(ITIs1<ITIMin)
    inx = ITIs1 > ITIMax  | ITIs1 < ITIMin;
    ITIs1(inx) = random('Normal',mng,sdg,1,sum(inx));
end
ITIs2 = random('Uniform',ITIMin,ITIMax,1,NumTrials);
prr = rand(1,NumTrials);
rr = zeros(2,NumTrials);
rr(1,prr<pmx1) = 1;
rr(2,prr>=pmx1) = 1;
ITIs = rr(1,:) .* ITIs1 + rr(2,:) .* ITIs2;
ITIMin = min(ITIs);
ITIMax = max(ITIs);

figure
hist(ITIs,50)

%% My exponential failure distribution

ITIMean = 1.4;
ITIMin = 0.1;
ITIMax = 5;
NumTrials = 10000;

ITIMean2 = ITIMean - ITIMin;
ITIMax2 = ITIMax - ITIMin;

temp = exprnd(ITIMean2,NumTrials,1);
while any(temp>ITIMax2)
    inx = temp > ITIMax2;
    temp(inx) = exprnd(ITIMean2,sum(inx),1)';
end
temp = temp + ITIMin;
ITIs = temp;

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
figure
for Phi = 0.26
intfg = zeros(size(cnts));
tau = cnts;
cntr = 1;   % counter
cnts2 = [cnts cnts(end)+dt:dt:cnts(end)+195*dt];
T = 0.3;
for t = cnts2;
    fg = ft .* exp(-(tau-t).^2/(2*(Phi*t^T)^2));
    intfg(cntr) = sum(fg*dt);
    cntr = cntr + 1;
end
fwt = 1 ./ (Phi * cnts2.^T * sqrt(2*pi)) .* intfg;
fwt = fwt / sum(fwt);

% Subjective hazard rate
Fwt = cumsum(fwt);
Rwt = 1 - Fwt;
hwt = (Rwt(1:end-1) - Rwt(2:end)) ./ (dt * Rwt(1:end-1));
hold on
plot(cnts2(1:end-1),hwt)
end