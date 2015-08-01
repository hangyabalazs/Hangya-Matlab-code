function simzshiftrun2
%SIMZSHIFTRUN   Runs SIMZSHIFT2.
%   SIMZSHIFTRUN calls SIMZSHIFT2 varying the LEN (segment length) parameter.
%
%   See also SIMZSHIFT2.

% Simulation
sno = 50;   % number of simulations
next = 1;
K = 5:5:30;
zm = zeros(1,length(K));
zmsd = zeros(1,length(K));
zmiqu = zeros(1,length(K));
zmiql = zeros(1,length(K));
zms = cell(1,length(K));
for k = K
    k
    lzm = zeros(1,sno);
    for t = 1:sno
        lzm(t) = simzshift2(1000,k,0.1,0.3,-50,4,20,'Poisson');
    end
    zms{next} = lzm;
    zm(next) = nanmedian(lzm);
    zmsd(next) = nanstd(lzm);
    zmiqu(next) = prctile(lzm(~isnan(lzm)),25);
    zmiql(next) = prctile(lzm(~isnan(lzm)),75);
    next = next + 1;
end

% Plot
figure
plot(K,zm,'g')
hold on
errorbar(K,zm,zmiql,zmiqu,'r+')
for k = 1:length(zms)
    lz = zms{k};
    lzn = lz(~isnan(lz));
    text(K(k)+0.005,zm(k)+5,num2str(length(lzn)))
end
1