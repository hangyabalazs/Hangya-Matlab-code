%% addition to FAlatjit.m for NB cells

tg = 6:15;   % from transgenic mice (ChAT-ChR2)
vrs=[1:5 16:22];   % from virus injected mice (ChAT-Cre)

%% punishment response

ranksum(FA_latency(tg),FA_latency(vrs))   % p = 0.4272
ranksum(FA_jitter(tg),FA_jitter(vrs))   % p = 0.3734
ranksum(FA_reliability(tg),FA_reliability(vrs))   % p = 0.3068

%% reward response: compare covariance controling for depth

depth = getvalue('DVpos',ChAT);
g = zeros(1,22);
g(vrs) = 1;
aoctool(depth,Hit_latency,g)  % see g*depth interaction: p = 0.65
aoctool(depth,Hit_jitter,g)  % see g*depth interaction: p = 0.12
aoctool(depth,Hit_reliability,g)   % see g*depth interaction: p = 0.64