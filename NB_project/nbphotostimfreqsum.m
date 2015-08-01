%% summary fig. for stim.-freq.-dependence of light response reliability

global DATAPATH
impdir = fullfile(DATAPATH,'NB\taggingsummary\reliability2\mat\');

dr = dir(impdir);
files = {dr(3:end).name};
% inx = [2 3 4 6 7 8 10 11 12 20 22 23 24 25];  indices for old file, when
% only behav sessions were included ('NB\taggingsummary\reliability\mat\')
inx = [2 3 5 7 8 9 11 12 13 23 25 26 28 29];
% 1 - 'n023_111218a_1.2' only 10 Hz stimulation
% 8 - 'n045_121217x_4.6' last point not good (drifting issue)
% 11 - 'n046_121230a_1.2' includes freely moving tagging session only to
% avoid drift effect!
inx = inx(3:end);   % first 2 cells: different set of frequencies
files = files(inx);
NumCells = length(files);   % n = 12
stim_freq = cell(1,NumCells);
Rs = nan(NumCells,4);
figure
hold on
for iC = 1:NumCells
    cellid = [files{iC}(1:14) '.' files{iC}(16)];
    load(fullfile(impdir,files{iC}))
    if cellid == 'n029_120214b_2.2'   % different frequencies
        Reliability(4) = NaN;
        burst_types(4) = 80;
    end
    figure
    plot(burst_types/2,Reliability)
    disp([num2str(iC) '  ' files{iC}])
%     keyboard
    stim_freq{iC} = burst_types / 2;
%     burst_types / 2
    Rs(iC,1:length(Reliability)) = Reliability;
end

% Mean
figure
plot([5 10 20 40],nanmean(Rs),'k','LineWidth',1)
hold on
E = errorbar([5 10 20 40],nanmean(Rs),nanse(Rs),'k','LineWidth',1);
errorbar_tick(E,0)   % eliminate horizontal line from errorbar
set(gca,'XLim',[4 41],'XTick',[5 10 20 40],'YLim',[0 0.7],'YTick',[0 0.7])
xlabel('Stimulation frequency')
ylabel('Spiking probablity')
setmyplot_balazs

%% summary fig. for stim.-freq.-dependence of light response reliability - with new cells

global DATAPATH
impdir = fullfile(DATAPATH,'NB\taggingsummary\reliability2\mat\');

dr = dir(impdir);
files = {dr(3:end).name};
% inx = [2 3 4 6 7 8 10 11 12 20 22 23 24 25];  indices for old file, when
% only behav sessions were included ('NB\taggingsummary\reliability\mat\')
inx = [2 3 5 7 8 9 11 12 13 23 25 26 28 29];
% 1 - 'n023_111218a_1.2' only 10 Hz stimulation
% 8 - 'n045_121217x_4.6' last point not good (drifting issue)
% 11 - 'n046_121230a_1.2' includes freely moving tagging session only to
% avoid drift effect!
inx = inx(3:end);   % first 2 cells: different set of frequencies
files = files(inx);
NumCells = length(files);   % n = 12
stim_freq = cell(1,NumCells);
Rs = nan(NumCells,4);
figure
hold on
for iC = 1:NumCells
    cellid = [files{iC}(1:14) '.' files{iC}(16)];
    load(fullfile(impdir,files{iC}))
    if cellid == 'n029_120214b_2.2'   % different frequencies
        Reliability(4) = NaN;
        burst_types(4) = 80;
    end
    figure
    plot(burst_types/2,Reliability)
    disp([num2str(iC) '  ' files{iC}])
%     keyboard
    stim_freq{iC} = burst_types / 2;
%     burst_types / 2
    Rs(iC,1:length(Reliability)) = Reliability;
end

impdir = fullfile(DATAPATH,'NB\taggingsummary_newdata\reliability2\mat\');
dr = dir(impdir);
files = {dr(3:end).name};
NumCells = length(files);   % n = 2
for iC = 1:NumCells
    cellid = [files{iC}(1:14) '.' files{iC}(16)];
    load(fullfile(impdir,files{iC}))
    figure
    plot(burst_types/2,Reliability)
    disp([num2str(iC) '  ' files{iC}])
%     keyboard
    stim_freq{end+1} = burst_types / 2;
%     burst_types / 2
    Rs(end+1,1:length(Reliability)) = Reliability;
end

impdir = fullfile(DATAPATH,'HDB\taggingsummary_newdata\reliability2\mat\');
dr = dir(impdir);
files = {dr(3:end).name};
NumCells = length(files);   % n = 8
for iC = 1:NumCells
    cellid = [files{iC}(1:14) '.' files{iC}(16)];
    load(fullfile(impdir,files{iC}))
    figure
    plot(burst_types/2,Reliability)
    disp([num2str(iC) '  ' files{iC}])
%     keyboard
    stim_freq{end+1} = burst_types / 2;
%     burst_types / 2
    Rs(end+1,1:length(Reliability)) = Reliability;
end

% Mean
figure
plot([5 10 20 40],nanmean(Rs),'k','LineWidth',1)
hold on
E = errorbar([5 10 20 40],nanmean(Rs),nanse(Rs),'k','LineWidth',1);
errorbar_tick(E,0)   % eliminate horizontal line from errorbar
set(gca,'XLim',[4 41],'XTick',[5 10 20 40],'YLim',[0 0.7],'YTick',[0 0.7])
xlabel('Stimulation frequency')
ylabel('Spiking probablity')
setmyplot_balazs

