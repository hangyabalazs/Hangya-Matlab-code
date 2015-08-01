%%

load('C:\Balazs\_data\VIP\freq_tuning_data_inh.mat')
load('C:\Balazs\_data\VIP\freq_tuning_data_VIP.mat')

%%

toneinx = 1:30;
NumTones = length(toneinx);
NumCells = length(inh);
[mu_stim, mu_nostim, max_stim, max_nostim] = deal(nan(1,NumCells));
for k = 1:NumCells;
    xc = toneinx;
    nostimrate = inh{k}.rate_tone;
    stimrate = inh{k}.rate_laser(:,k);
%     naninx = isnan(stimrate) | isnan(nostimrate);
%     nostimrate(naninx) = [];
%     stimrate(naninx) = [];
%     xc(naninx) = [];
    nninx = find(~isnan(stimrate)&~isnan(nostimrate));
    nostimrate = linterp(xc(nninx),nostimrate(nninx),xc);
    stimrate = linterp(xc(nninx),stimrate(nninx),xc);
    figure
    plot(xc,nostimrate,'ko')
    hold on
    plot(xc,stimrate,'o','MarkerFaceColor','b')

    mloc1 = find(nostimrate==max(nostimrate));
    mloc1 = mloc1(1);
    mloc2 = find(stimrate==max(stimrate));
    mloc2 = mloc2(1);
    mloc = mean(mloc1,mloc2);
    
    sstim = smooth(stimrate,'linear',3);
    snostim = smooth(nostimrate,'linear',3);
    
    s = fitoptions('Method','NonlinearLeastSquares',...
        'Lower',[0,0,0,0],...
        'Upper',[Inf,Inf,xc(end),Inf],...
        'Startpoint',[max(snostim),0,xc(mloc) 1],...
        'Robust','on');
    f = fittype('a*(sqrt(2*pi*sigma^2)^(-1)*exp(-(x-mu)^2/(2*sigma^2)))+b','options',s);
    [fun gof] = fit(xc',snostim,f,s);
    mu_nostim(k) = fun.mu;
    
    hold on
    plot(xc,snostim,'k:')
    L = plot(fun,'k');
    max_nostim(k) = max(get(L,'YData'));
    
    s = fitoptions('Method','NonlinearLeastSquares',...
        'Lower',[0,0,0,0],...
        'Upper',[Inf,Inf,xc(end),Inf],...
        'Startpoint',[max(sstim),0,xc(mloc) 1],...
        'Robust','on');
    f = fittype('a*(sqrt(2*pi*sigma^2)^(-1)*exp(-(x-mu)^2/(2*sigma^2)))+b','options',s);
    [fun gof] = fit(xc',sstim,f,s);
    mu_stim(k) = fun.mu;
    
    hold on
    plot(xc,sstim,'b:')
    L = plot(fun,'b');
    max_stim(k) = max(get(L,'YData'));

end

%%

[stim_aligned, nostim_aligned] = deal(nan(NumCells,NumTones*2));
for k = 1:NumCells
    nostimrate = inh{k}.rate_tone;
    stimrate = inh{k}.rate_laser(:,k);
    nninx = find(~isnan(stimrate)&~isnan(nostimrate));
    nostimrate = linterp(xc(nninx),nostimrate(nninx),xc);
    stimrate = linterp(xc(nninx),stimrate(nninx),xc);
    
    lside = sum(xc<mu_stim(k));
    rside = NumTones - lside;
    stim_aligned(k,(NumTones-lside):(NumTones+rside-1)) = stimrate;
    
    lside = sum(xc<mu_nostim(k));
    rside = NumTones - lside;
    nostim_aligned(k,(NumTones-lside):(NumTones+rside-1)) = nostimrate;
end

%%

figure
plot(nanmean(nostim_aligned),'k')
hold on
plot(nanmean(stim_aligned))
errorshade(1:2*NumTones,nanmean(nostim_aligned),nanstd(nostim_aligned)./sqrt(sum(~isnan(nostim_aligned))),...
    'ShadeColor','k','LineColor','k')
errorshade(1:2*NumTones,nanmean(stim_aligned),nanstd(stim_aligned)./sqrt(sum(~isnan(stim_aligned))),...
    'ShadeColor','k','LineColor','b')

%%

mx = max(nostim_aligned,[],2);

figure
plot(nanmean(nostim_aligned./repmat(mx,1,size(nostim_aligned,2))),'k')
hold on
plot(nanmean(stim_aligned./repmat(mx,1,size(nostim_aligned,2))))
% errorshade(1:2*NumTones,nanmean(nostim_aligned),nanstd(nostim_aligned)./sqrt(sum(~isnan(nostim_aligned))),...
%     'ShadeColor','k','LineColor','k')
% errorshade(1:2*NumTones,nanmean(stim_aligned),nanstd(stim_aligned)./sqrt(sum(~isnan(stim_aligned))),...
%     'ShadeColor','k','LineColor','b')

%%

mx = max(nostim_aligned,[],2);   % norm. to max. response
nostim_aligned_norm = nostim_aligned ./ repmat(mx,1,size(nostim_aligned,2));
stim_aligned_norm = stim_aligned ./ repmat(mx,1,size(stim_aligned,2));
nostim_mn = nanmean(nostim_aligned_norm);
stim_mn = nanmean(stim_aligned_norm);
nostim_err = nanstd(nostim_aligned_norm)./sqrt(sum(~isnan(nostim_aligned_norm)));
stim_err = nanstd(stim_aligned_norm)./sqrt(sum(~isnan(stim_aligned_norm)));
xc = 1:length(nostim_mn);
naninx = isnan(nostim_mn) | isnan(stim_mn);
nostim_mn(naninx) = [];
stim_mn(naninx) = [];
nostim_err(naninx) = [];
stim_err(naninx) = [];
xc(naninx) = [];

figure
% plot(xc,nostim_mn,'ko')
errorbar(xc,nostim_mn,nostim_err,'ko')
hold on
% plot(xc,stim_mn,'bo')
errorbar(xc,stim_mn,stim_err,'bo')

mloc1 = find(nostim_mn==max(nostim_mn));
mloc1 = mloc1(1);
mloc2 = find(stim_mn==max(stim_mn));
mloc2 = mloc2(1);
mloc = mean(mloc1,mloc2);

s = fitoptions('Method','NonlinearLeastSquares',...
    'Lower',[0,0,0,0],...
    'Upper',[Inf,Inf,xc(end),Inf],...
    'Startpoint',[max(stim_mn),0,xc(mloc) 1],...
    'Robust','on');
f = fittype('a*(sqrt(2*pi*sigma^2)^(-1)*exp(-(x-mu)^2/(2*sigma^2)))+b','options',s);
[fun gof] = fit(xc',stim_mn',f,s);

hold on
plot(fun,'b')

s = fitoptions('Method','NonlinearLeastSquares',...
    'Lower',[0,0,0,0],...
    'Upper',[Inf,Inf,xc(end),Inf],...
    'Startpoint',[max(nostim_mn),0,xc(mloc) 1],...
    'Robust','on');
f = fittype('a*(sqrt(2*pi*sigma^2)^(-1)*exp(-(x-mu)^2/(2*sigma^2)))+b','options',s);
[fun gof] = fit(xc',nostim_mn',f,s);

hold on
plot(fun,'k')

%%

% mx = nostim_aligned(:,NumTones);   % Scanziani (norm. to peak response)
mx = max_nostim';   % norm to peak of the fit
nostim_aligned_norm = nostim_aligned ./ repmat(mx,1,size(nostim_aligned,2));
stim_aligned_norm = stim_aligned ./ repmat(mx,1,size(stim_aligned,2));
nostim_mn = nanmean(nostim_aligned_norm);
stim_mn = nanmean(stim_aligned_norm);
nostim_err = nanstd(nostim_aligned_norm)./sqrt(sum(~isnan(nostim_aligned_norm)));
stim_err = nanstd(stim_aligned_norm)./sqrt(sum(~isnan(stim_aligned_norm)));
xc = 1:length(nostim_mn);
naninx = isnan(nostim_mn) | isnan(stim_mn);
nostim_mn(naninx) = [];
stim_mn(naninx) = [];
nostim_err(naninx) = [];
stim_err(naninx) = [];
xc(naninx) = [];

figure
% plot(xc,nostim_mn,'ko')
errorbar(xc,nostim_mn,nostim_err,'ko')
hold on
% plot(xc,stim_mn,'bo')
errorbar(xc,stim_mn,stim_err,'bo')

mloc1 = find(nostim_mn==max(nostim_mn));
mloc1 = mloc1(1);
mloc2 = find(stim_mn==max(stim_mn));
mloc2 = mloc2(1);
mloc = mean(mloc1,mloc2);

s = fitoptions('Method','NonlinearLeastSquares',...
    'Lower',[0,0,0,0],...
    'Upper',[Inf,Inf,xc(end),Inf],...
    'Startpoint',[max(stim_mn),0,xc(mloc) 1],...
    'Robust','on');
f = fittype('a*(sqrt(2*pi*sigma^2)^(-1)*exp(-(x-mu)^2/(2*sigma^2)))+b','options',s);
[fun gof] = fit(xc',stim_mn',f,s);

hold on
plot(fun,'b')

s = fitoptions('Method','NonlinearLeastSquares',...
    'Lower',[0,0,0,0],...
    'Upper',[Inf,Inf,xc(end),Inf],...
    'Startpoint',[max(nostim_mn),0,xc(mloc) 1],...
    'Robust','on');
f = fittype('a*(sqrt(2*pi*sigma^2)^(-1)*exp(-(x-mu)^2/(2*sigma^2)))+b','options',s);
[fun gof] = fit(xc',nostim_mn',f,s);

hold on
plot(fun,'k')