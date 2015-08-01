%% parameters

dt = 0.01;  % delta t
t0 = 0;   % time 0
te = 10;  % last evaluated time point
t = t0:dt:te;   % time vector
tk = [1 1.5 6];   % event times (members of the time vector)
lambda = 1.5;   % time constant
N0 = 2;   % step size at event times

%% simulation

tic
cntr = 1;
L1t = nan(1,length(t)-1);   % simulation in time steps
L2t = nan(1,length(t)-1);   % simulation using the closed form
L1t(1) = 0;   % t=0
L2t(1) = 0;   % t=0
for ti = t(2:end)
    cntr = cntr + 1;
    if ~ismember(ti,tk);
        L1t(cntr) = L1t(cntr-1) * exp(-lambda*dt);  % exponential decay
    else
        L1t(cntr) = L1t(cntr-1) * exp(-lambda*dt) + N0;   % an event occurred
    end
    
    tki = tk(tk<=ti);   % events before ti
    L2t(cntr) = N0 * sum(exp(-lambda*(ti-tki)));
end
toc

%% evaluation

L2te = N0 * sum(exp(-lambda*(te-tk)));   % closed form, last time point
L2te - L1t(end)

figure
plot(t,L1t)
hold on
plot(t,L2t,'r')

%% calculate with convolution

tic
z = zeros(size(t));
z(round(tk/dt)+1) = N0;
flt = exp(-lambda*(t));
L3t = conv2(z,flt);
L3t = L3t(1:length(t));
toc

figure
plot(t,L1t)
hold on
plot(t,L3t,'g')