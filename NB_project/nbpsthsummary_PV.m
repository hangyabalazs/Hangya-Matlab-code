function nbpsthsummary_PV
%NBPSTHSUMMARY_PV   PSTH summary figures.
%   NBPSTHSUMMARY_PV plots population PSTH for activated PV+ cells.
%   Z-scored PSTHs are used throughout the analysis.
%
%   Average PSTH aligned to reward or punishment is plotted for all PV+
%   cells. A heat map of the population PSTH for is also plotted separately
%   for hit and false alarm trials.
%
%   See also NBPSTHSUMMARY.

%   Balazs Hangya, Cold Spring Harbor Laboratory
%   1 Bungtown Road, Cold Spring Harbor
%   balazs.cshl@gmail.com
%   2-Oct-2013

%   Edit log: BH 10/2/13

% Directories
global DATAPATH
resdir = ([DATAPATH 'NB\psthsummary_PV\']);   % result directory

% Groups of identified NB cells
selstr = ['"PV+"==1&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})'];
PV = selectcell(selstr);   % cell IDs for PV cells

% Average PSTHs for identified PV cells
green = [51 204 51] / 255;   % colors for plotting
red = [216 41 0] / 255;
NumPV = length(PV);   % number of identified PV+ neurons
FA_psth = getvalue('FA_psth',PV);    % PSTH aligned to air puff 
FA_psth = nancell2mat(FA_psth);
Hit_psth = getvalue('Hit_psth',PV);   % PSTH aligned to water delivery
Hit_psth = nancell2mat(Hit_psth);
for k = 1:NumPV   % z-score PSTHs
    FA_psth(k,:) = (FA_psth(k,:) - mean(Hit_psth(k,:))) / std(Hit_psth(k,:));
    Hit_psth(k,:) = zscore(Hit_psth(k,:));
end

figure   % plot
time = -600:600;
FA_psth = FA_psth(~all(isnan(FA_psth),2),:);
Hit_psth = Hit_psth(~all(isnan(Hit_psth),2),:);
errorshade(time,mean(FA_psth),std(FA_psth)/sqrt(size(FA_psth,1)),...
    'LineColor',red,'ShadeColor',red)
hold on
errorshade(time,mean(Hit_psth),std(Hit_psth)/sqrt(size(Hit_psth,1)),...
    'LineColor',green,'ShadeColor',green)
saveas(gcf,fullfile(resdir,'PV_psth_Hit_FA_average.fig'))

% Population PSTH
figure   % plot all FA PSTHs
imagesc(time,1:size(FA_psth,1),FA_psth)
colormap(hot)
saveas(gcf,fullfile(resdir,'PV_poppsth_FA_all.fig'))

[srt Ia] = sort(max(Hit_psth,[],2),'descend');   % sort based on Hit response
figure   % plot all FA PSTHs, sorted
imagesc(time,1:size(FA_psth,1),FA_psth(Ia,:))
colormap(hot)
saveas(gcf,fullfile(resdir,'PV_poppsth_FA_all_sorted.fig'))

figure   % plot all Hit PSTHs
imagesc(time,1:size(Hit_psth,1),Hit_psth)
colormap(hot)
saveas(gcf,fullfile(resdir,'PV_poppsth_Hit_all.fig'))

figure   % plot all Hit PSTHs; sort based on Hit response
imagesc(time,1:size(Hit_psth,1),Hit_psth(Ia,:))
colormap(hot)
saveas(gcf,fullfile(resdir,'PV_poppsth_Hit_all_sorted.fig'))

keyboard