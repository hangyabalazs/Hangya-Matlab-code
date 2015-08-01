function simzshiftrun4
%SIMZSHIFTRUN   Runs SIMZSHIFT2.
%   SIMZSHIFTRUN calls SIMZSHIFT2 and plots Z-shift distribution.
%
%   See also SIMZSHIFT2.

% Simulation
sno = 1000;     % number of simulations
lzm = zeros(1,sno);
for t = 1:sno
    t
    lzm(t) = simzshift2(1000,10,0.1,0.3,-50,4,20,'Poisson');
end

% Plot
figure
hist(lzm,50)
1