%% phasediff vs. rt
% Directories
global DATADIR
global DATAPATH
inpdir = [DATADIR 'human_SG\'];
resdir = [DATAPATH 'Expectancy\phasediff\'];
cd(resdir)

% Import
ff = [inpdir 'delta_theta_phasediff_rt.mat'];
load(ff)
data = data_phasediff_rt;
data(data<-pi) = data(data<-pi) + 2 * pi;
data(data>pi) = data(data>pi) - 2 * pi;

dim2 = 4;       % condition: 91%
dim3 = 1;       % channel: Fz
phasediff = cell(1,13);
rt = cell(1,13);
R = cell(1,13);
for k = 1:13
    phasediff{k} = squeeze(data(k,dim2,dim3,1,:));
    rt{k} = squeeze(data(k,dim2,dim3,2,:));
    [b,bint,r,rint,stats] = regress(rt{k},...
        [ones(length(phasediff{k}),1) sin(phasediff{k}) cos(phasediff{k})]);
    R{k} = sqrt(stats(1));
    H = figure;
    plot(phasediff{k},rt{k},'.')
    x_lim = xlim;
    y_lim = ylim;
    text(x_lim(1)+(x_lim(2)-x_lim(1))*0.3,y_lim(1)+(y_lim(2)-y_lim(1))*0.9,...
        ['circular-linear correlation: ' num2str(R{k})],'Color','red')
    saveas(H,['phasediff' num2str(k) '.fig'])
end

close all

%% phasediff vs. rt 2

cd(resdir)

dim2 = 4;       % condition: 91%
dim3 = 1;       % channel: Fz
for k = 1:13
    phasediff_vs_rt = squeeze(data(k,dim2,dim3,1:2,:))';
    spvr = sortrows(phasediff_vs_rt,1);
    H = figure;
    subplot(2,1,1)
    plot(spvr(:,1))
    subplot(2,1,2)
    plot(spvr(:,2))
    saveas(H,['sorted_phasediff' num2str(k) '.fig'])
end
close all

%% phasediff hist.

cd(resdir)

edges = -pi:2*pi/9:pi;
cnts = (edges(1:end-1) + edges(2:end)) / 2;
phasediff = cell(1,13);
for k = 1:13
    phasediff{k} = squeeze(data(k,dim2,dim3,1,:));
    [nm xout] = histc(phasediff{k},edges);
    nm = nm(1:end-1);
    H = figure;
    plot([cnts cnts+2*pi],[nm' nm'],'r')
    saveas(H,['phasediff_hist' num2str(k) '.fig'])
end
close all