%% parameters

EL = -65;  % in mV
Rm = 10;  % in MOhms
tau = 10;   % in ms
Vthr = -50;   % mV
Vreset = -65;  % mV

%% time

dt = 0.1;
Tend = 500;
T = 0+dt:dt:Tend;

%% Ie

z = zeros(1,Tend/dt);
% ps = round(rand(1,10)*Tend/dt);
ps = [50 100 180 220 300 350 380 450] / dt;
lenp = length(ps);
pri = nan(lenp,Tend/dt);
for k = 1:lenp
    z2 = z;
    z2(ps(k)) = 1;
    gwn = gausswin(501,rand*9+1);
    gwn = gwn / sum(gwn);
    pri(k,:) = filtfilt(gwn,1,z2);
end
Ie = sum(pri) * 500;

%% V

Vinf = EL + Rm * Ie;
V0 = EL;

V = nan(1,Tend/dt);
V(1) = Vinf(1) + (V0 - Vinf(1)) * exp(-dt/tau);
for k = 2:Tend/dt
    if V(k-1) > Vthr;
        V(k-1) = 0;
        V(k) = Vreset;
        continue
    end
    V(k) = Vinf(k) + (V(k-1) - Vinf(k)) * exp(-dt/tau);
end