%%

d_hat = 0:0.001:15;
mu = 0;
sigma = 1;
I = trapz(d_hat,normpdf(d_hat,mu,sigma))   % normalize the percept distribution to integral 1
p_d_hat = 1 / I * normpdf(d_hat,mu,sigma);
c_d_hat = normcdf(d_hat,mu,sigma);
trapz(d_hat,p_d_hat.*c_d_hat)

%%
figure
plot(d_hat,p_d_hat.*c_d_hat)
figure;plot(d_hat,c_d_hat)
figure;plot(d_hat,p_d_hat)
figure;plot(c_d_hat,p_d_hat,'.k')   % confidence-probability pairs

%%

dst = 'Normal';   % noise distribution
d_hat = 0:0.001:15;
mu = 0;
sigma = 1;
I = trapz(d_hat,pdf(dst,d_hat,mu,sigma))   % normalize the percept distribution to integral 1
p_d_hat = 1 / I * pdf(dst,d_hat,mu,sigma);
c_d_hat = cdf(dst,d_hat,mu,sigma);
trapz(d_hat,p_d_hat.*c_d_hat)

%%

dst = 't';   % noise distribution
d_hat = 0:0.001:15;
mu = 15;
sigma = 1;
I = trapz(d_hat,pdf(dst,d_hat,mu,sigma))   % normalize the percept distribution to integral 1
p_d_hat = 1 / I * pdf(dst,d_hat,mu,sigma);
c_d_hat = cdf(dst,d_hat,mu,sigma);
trapz(d_hat,p_d_hat.*c_d_hat)

%%

dst = 'Normal';   % noise distribution
d_hat = 0:0.001:15;
mu = 0;
sigma = 1;
I = trapz(d_hat,pdf(dst,d_hat,mu,sigma))   % normalize the percept distribution to integral 1
p_d_hat = 1 / I * pdf(dst,d_hat,mu,sigma);
c_d_hat = cdf(dst,d_hat,mu,sigma);
trapz(d_hat,p_d_hat.*c_d_hat)

%%

dst = 'Normal';   % noise distribution
d_hat = -15:0.001:0;
mu = 0;
sigma = 2;
trapz(d_hat,cdf(dst,d_hat,mu,sigma)) 