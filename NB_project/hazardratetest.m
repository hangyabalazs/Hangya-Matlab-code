%% Failure distribution

ITIMean = 4;
ITIMin = 3;
ITIMax = 9;
NumTrials = 10000;

temp = exprnd(ITIMean,NumTrials,1);
temp(temp<ITIMin) = ITIMin;
temp(temp>ITIMax) = ITIMax;
ITIs = temp;

%% New failure distribution

ITIMean = 1.5;
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

dt = 0.1;
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
intfg = zeros(size(cnts));
tau = cnts;
cntr = 1;   % counter
cnts2 = [cnts(1)-95*dt:dt:cnts(1)-dt cnts cnts(end)+dt:dt:cnts(end)+195*dt];
for t = cnts2
    fg = ft' .* exp(-(tau-t).^2/(2*Phi^2*t^2));
    intfg(cntr) = sum(fg*dt);
    cntr = cntr + 1;
end
fwt = 1 ./ (Phi * cnts2 * sqrt(2*pi)) .* intfg;
sum(fwt)
fwt = fwt / sum(fwt);

% Subjective hazard rate
Fwt = cumsum(fwt);
Rwt = 1 - Fwt;
hwt = (Rwt(1:end-1) - Rwt(2:end)) ./ (dt * Rwt(1:end-1));
figure
plot(cnts2(1:end-1),hwt)