function L5spiking5percent_call
%L5SPIKING5PERCENT_CALL    Calls L5_SPIKING_5PERCENT for a sequence of directories.
%   L5SPIKING5PERCENT_CALL calculates histogram of pooled cycle first spike
%   data and corresponding 5 and 10 percentile values for layer 5 pyramidal
%   cells.
%
%   Note: slow oscillatoin cut at 500 ms (2 Hz)!
%
%   See also L5_SPIKING_5PERCENT.

% Directories
global DATAPATH
inproot = [DATAPATH 'Hajni_layer_5\disc2\'];
dr = dir(inproot);
inpadd = {};
for k = 3:length(dr)
    if dr(k).isdir
        inpadd{end+1} = dr(k).name;
    end
end

% Call L5_SPIKING_5PERCENT
aang = [];
aang_cfsp = [];
for k = 1:length(inpadd)
    inpdir = [inproot inpadd{k} '\']
    [paang paang_cfsp resdir1] = L5_spiking_5percent(inpdir);
    aang = [aang paang];
    aang_cfsp = [aang_cfsp paang_cfsp];
end

% Phase histogram for cycle first spike lower 10%
edges = -180:20:180;     % edges for phase histogram
cnts = (edges(1:end-1) + edges(2:end)) / 2;
n_cfsp = length(aang_cfsp);
ftm_cfsp = sum(exp(1).^(i*aang_cfsp)) / n_cfsp;    % first trigonometric moment
ang_cfsp = angle(ftm_cfsp);   % mean angle
mvl_cfsp = abs(ftm_cfsp);     % mean resultant length
aang_cfsp = aang_cfsp * 180 / pi;
ang_cfsp = ang_cfsp * 180 / pi;
[nm_cfsp,xout_cfsp] = histc(aang_cfsp,edges);   % phase histogram
nm_cfsp = nm_cfsp(1:end-1);

saang_cfsp = sort(aang_cfsp,'ascend');
l5 = length(aang_cfsp) * 0.10;
aang_cfsp5 = saang_cfsp(1:round(l5));
n_cfsp5 = length(aang_cfsp5);
% ftm_cfsp5 = sum(exp(1).^(i*aang_cfsp5/180*pi)) / n_cfsp5;    % first trigonometric moment
% ang_cfsp5 = angle(ftm_cfsp5);   % mean angle
% mvl_cfsp5 = abs(ftm_cfsp5);     % mean resultant length
% aang_cfsp5 = aang_cfsp5 * 180 / pi;
% ang_cfsp5 = ang_cfsp5 * 180 / pi;
[nm_cfsp5,xout_cfsp5] = histc(aang_cfsp5,edges);   % phase histogram
nm_cfsp5 = nm_cfsp5(1:end-1);

% Display percentiles
pct5 = prctile(aang_cfsp,5);
disp(['5%: ' num2str(pct5)])
pct10 = prctile(aang_cfsp,10);
disp(['10%: ' num2str(pct10)])

% Plot
H1 = figure;      % cycle first spikes
bar(cnts,nm_cfsp5'/n_cfsp5)
title(gca,'Cycle first spike for 5 L5 cells, first 10 percents')
y_lim = ylim;
axis([-200 200 y_lim(1) y_lim(2)])
% str = ['\it{Mean angle: }' '\bf ' num2str(ang_cfsp5)];
% text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','red')
% str = ['\it{Mean resultant length: }' '\bf ' num2str(mvl_cfsp5)];
% text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/12,str,'Color','red')
% str = ['\it{n: }' '\bf ' num2str(n_cfsp5)];
% text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/6,str,'Color','red')
set(gca,'XTick',[-180 -90 0 90 180])

H2 = figure;
bar(cnts,nm_cfsp'/n_cfsp)
title(gca,['Cycle first spike for 5 L5 cells'])
y_lim = ylim;
axis([-200 200 y_lim(1) y_lim(2)])
str = ['\it{Mean angle: }' '\bf ' num2str(ang_cfsp)];
text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','red')
str = ['\it{Mean resultant length: }' '\bf ' num2str(mvl_cfsp)];
text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/12,str,'Color','red')
str = ['\it{n: }' '\bf ' num2str(n_cfsp)];
text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/6,str,'Color','red')
set(gca,'XTick',[-180 -90 0 90 180])

H3 = figure;
edges = -180:10:180;     % edges for phase histogram
cnts = (edges(1:end-1) + edges(2:end)) / 2;
[nm_cfsp,xout_cfsp] = histc(aang_cfsp,edges);   % phase histogram
nm_cfsp = nm_cfsp(1:end-1);plot(cnts,nm_cfsp'/n_cfsp,'k')
hold on
ml5 = prctile(aang_cfsp,10);
xml5 = linterp(cnts,nm_cfsp'/n_cfsp,ml5);
area([cnts(cnts<ml5) ml5],[nm_cfsp(cnts<ml5)/n_cfsp xml5],'FaceColor','red')
title(gca,['Cycle first spike for 5 L5 cells'])
y_lim = ylim;
axis([-200 200 y_lim(1) y_lim(2)])
str = ['\it{Mean angle: }' '\bf ' num2str(ang_cfsp)];
text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6,str,'Color','red')
str = ['\it{Mean resultant length: }' '\bf ' num2str(mvl_cfsp)];
text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/12,str,'Color','red')
str = ['\it{n: }' '\bf ' num2str(n_cfsp)];
text(60,y_lim(2)-(y_lim(2)-y_lim(1))/6-(y_lim(2)-y_lim(1))/6,str,'Color','red')
set(gca,'XTick',[-180 -90 0 90 180])

% Save
cd(resdir1)
fns = 'CYCLEFIRST10.fig';
saveas(H1,fns)
fns = 'CYCLEFIRST10.jpg';
saveas(H1,fns)
fns = 'CYCLEFIRSTall.fig';
saveas(H2,fns)
fns = 'CYCLEFIRSTall.jpg';
saveas(H2,fns)
fns = 'CYCLEFIRSTall_10.fig';
saveas(H3,fns)
fns = 'CYCLEFIRSTall_10.jpg';
saveas(H3,fns)
save 'CYCLEFIRST.mat' 'aang' 'aang_cfsp' 'pct5' 'pct10'
close all