%% My bimodal failure distribution

NumTrials = 10000;
ITIMin = 0.1;
ITIMax = 3;
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

%% My bimodal failure distribution

NumTrials = 10000;
ITIMin = -10;
ITIMax = 10;
mng1 = 0.3;   % parameters for the Gaussians
mng2 = 2;
sdg = 0.15;
pmx1 = 0.5;   % mixing probabilities
pmx2 = 0.5;
pmx3 = 1 - pmx1 - pmx2;

ITIs1 = random('Normal',mng1,sdg,1,NumTrials);
ITIs2 = random('Normal',mng2,sdg,1,NumTrials);
ITIs3 = random('Uniform',0.1,3,1,NumTrials);
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

%%


[m,s,Pa,iter,Q_tot,e_tot] = em_alg_function(ITIs,[0.3 2],[0.15^2 0.15^2],[0.5 0.5],0.1);

%%

figure;
[nm xout] = hist(ITIs,100);
bar(xout,nm/sum(nm))

gs = Pa(1) * normpdf(xout,m(1),sqrt(s(1))) + Pa(2) * normpdf(xout,m(2),sqrt(s(2)));
hold on
plot(xout,gs/sum(gs),'r')