% Probabilities
% p = [2 3 2 6 6 7 6 8 10 10 9 12 11 12 13 14 15 14 13 12 11 10 10 ...
%     9 8 7 7 6 6 5 5 5 5 5 5 5 5 5 5 5 5 5] / 100;
st = 500;   % stim. time
p0 = ones(1,st) * 0.01;     % baseline prob.
rp = [0.15 0.20 0.05];   % response prob.
p = [p0 ones(1,3)*0.01 rp p0 p0] / 5;
realrp = rp(1);
for k = 2:length(rp)
    realrp = realrp +(1 - realrp) * rp(k);  % real (global) response prob.
end
tl = length(p);

% Simulate trials
tno = 100;   % number of trials
spt = zeros(tno,tl);
for k = 1:tno
    sptr = double(rand(1,tl)<p);    % sample independendtly
    
    ds = sptr(1:end-1) .* sptr(2:end);   % eliminate spikes within <3 distance
    fds = find(ds,1,'first');
    while ~isempty(fds)
        sptr(fds+1) = 0;
        ds = sptr(1:end-1) .* sptr(2:end);
        fds = find(ds,1,'first');
    end
    ds = sptr(1:end-2) .* sptr(3:end);
    fds = find(ds,1,'first');
    while ~isempty(fds)
        sptr(fds+2) = 0;
        ds = sptr(1:end-2) .* sptr(3:end);
        fds = find(ds,1,'first');
    end
    spt(k,:) = sptr;
end

% Raster plot
figure
imagesc(spt)
line([st st],[0 tl],'Color','red')
C = colormap('bone');
colormap(flipud(C))

% Probability fit
figure
plot(p,'r')
hold on
plot(sum(spt./tno),'k')

% Probability ratio
wn = 10;    % window size
pa = sum(sum(spt(:,st:st+wn))) ./ (tno*wn);   % firing prob. after stim.
pb = sum(sum(spt(:,st-wn:st))) ./ (tno*wn);   % firing prob. before stim.
prob_ratio =  pa / pb;

% Stimulus 'postdiction'
kis = zeros(tno,tl);
jt = zeros(1,tl);
jt2 = zeros(1,tl);
lt = zeros(1,tl);
lt2 = zeros(1,tl);
for t = 1:tl
    for k = 1:tno
        cspt = spt(k,t:end);
        pki = find(cspt,1,'first');
        if ~isempty(pki)
            kis(k,t) = pki;
        else
            kis(k,t) = NaN;
        end
    end
    jt(t) = nanstd(kis(:,t));   % jitter
    nkis = kis(:,t);
    nkis = nkis(~isnan(nkis));
    jt2(t) = iqr(nkis);
    lt(t) = nanmean(kis(:,t));   % latency
    lt2(t) = nanmedian(kis(:,t));
end
figure
plot(jt(1:1000))
figure
plot(jt2(1:1000))
figure
plot(lt(1:1000))
figure
plot(lt2(1:1000))

% KL-distance
nm = floor((tl-1)/10);
lsi = zeros(tno,nm);
slsi = zeros(tno,nm);
hlsi = zeros(11,nm);
nhlsi = zeros(11,nm);   % normalized ISI histogram
next = 1;
for t = 1:10:tl-9
    for k = 1:tno
        cspt = spt(k,t:t+9);
        pki = find(cspt,1,'first');
        if ~isempty(pki)
            lsi(k,next) = pki;
        else
            lsi(k,next) = 0;
        end
    end
    slsi(:,next) = sort(lsi(:,next));
    hlsi(:,next) = hist(slsi(:,next),11);
    nhlsi(:,next) = hlsi(:,next) / sum(hlsi(:,next));
    next = next + 1;
end
figure      % plot ISIs
imagesc(lsi)
figure      % plot sorted ISIs
imagesc(slsi)

kn = st / 10 + 1;   % pairwise KL-distance
kld = zeros(kn,kn);
for k1 = 1:kn
    for k2 = 1:kn
        kld(k1,k2) = KLdist(nhlsi(:,k1),nhlsi(:,k2));
    end
end
figure    % plot KL-distance
imagesc(kld)