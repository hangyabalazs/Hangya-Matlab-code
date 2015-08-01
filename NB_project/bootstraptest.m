%% Normal sample

sample_dist = randn(1,100);
n = numel(sample_dist);
mn = mean(sample_dist);
se = std(sample_dist) / sqrt(n);
se_median = se * 1.253

%% Bootstrap SE

bno = 10000;
bootstrap_mn = nan(1,bno);
bootstrap_mn2 = nan(1,bno);
for k = 1:bno
    bootstrap_sample = sample_dist(randi(n,1,n));
    bootstrap_mn(k) = mean(bootstrap_sample);
    bootstrap_mn2(k) = median(bootstrap_sample);
end
bootstrap_se = std(bootstrap_mn);
bootstrap_se_median = std(bootstrap_mn2)