%% load

% Import session data from Excel
global DATAPATH
tblfile = [DATAPATH 'SniffWhisk\data_selection.xlsx'];
[tbl0 tbl] = xlsread(tblfile);

% Call 'snwdisc'
numSessions = size(tbl,1);
allsegments = [];
sr = 1000;
for iS = 1:numSessions
    rrat = tbl{iS,1};
    rdate = tbl{iS,2};
    wl = tbl0(iS,3);
    
    % Load data for phase calculation
    fnme = [DATAPATH 'SniffWhisk\phasehist8\PHASE_' rrat '_' rdate '.mat'];
    load(fnme)
    
    for k = 1:length(segments)
        segments(k).start = segments(k).start * wl / sr;
        segments(k).end = segments(k).end * wl / sr;
        segments(k).inhalation_start = segments(k).inhalation_start * wl / sr;
        segments(k).inhalation_end = segments(k).inhalation_end * wl / sr;
        segments(k).exhalation_start = segments(k).exhalation_start * wl / sr;
        segments(k).exhalation_end = segments(k).exhalation_end * wl / sr;
        segments(k).whisk_events = segments(k).whisk_events * wl / sr;
        segments(k).raw_whisk_events = segments(k).raw_whisk_events * wl / sr;
    end
    
    allsegments = [allsegments segments];
end
segments = allsegments;


%% filter

segments2 = segments(1);
for k = 1:length(segments)
    if segments(k).length > 1 && segments(k).group == 1
        segments2(end+1) = segments(k);
    end
end
segments2 = segments2(2:end);

%% sort

lense = length(segments2);
fr = nan(1,lense);
for k = 1:lense
   fr(k) = segments2(k).frequency;
end
[sr ins] = sort(fr);

%% add more

fld = fieldnames(segments);
tmpp = segments2(235);
for k = 6:10
    tmp = tmpp.(fld{k});
    tmpp.(fld{k}) = tmp(1:end);
end
tmpp.start = tmpp.inhalation_start(1);
segments3 = [segments3(1:5) tmpp segments3(6:end)];


%% plot

figure
cntr = 0;

% inss = ins([1 4 10 13]);
% inss = ins;

T = setdiff(1:length(segments3),[132 174 234 530 613]);
% T = 1:36; % checked examples to 100


for k = T
    try
    cntr = cntr + 1;
    ln = length(segments3(k).inhalation_start);
    figpos = get(gcf,'Position');
    X = [segments3(k).inhalation_start' segments3(k).inhalation_end' ...
        segments3(k).inhalation_end' segments3(k).inhalation_start']-segments3(k).start;
    Y = repmat([cntr cntr cntr+1 cntr+1],ln,1);
    patch(X',Y',[1 1 0.5],'EdgeColor','none')
    
    line([segments3(k).whisk_events; segments3(k).whisk_events]-segments3(k).start,...
        [ones(1,length(segments3(k).whisk_events))*cntr; ...
        ones(1,length(segments3(k).whisk_events))*cntr+1],...
        'LineWidth',3,'Color',[0 153 0]/255)
    end
end

% xlim([0 55])

%% plot seg4

figure
cntr = 0;

% inss = ins([1 4 10 13]);
% inss = ins;

T = 500:550;

for k = T
    cntr = cntr + 1;
    ln = length(segments4(k).inhalation_start);
    figpos = get(gcf,'Position');
    X = [segments4(k).inhalation_start' segments4(k).inhalation_end' ...
        segments4(k).inhalation_end' segments4(k).inhalation_start']-segments4(k).start;
    Y = repmat([cntr cntr cntr+1 cntr+1],ln,1);
    patch(X',Y',[1 1 0.5])
    
    line([segments4(k).whisk_events; segments4(k).whisk_events]-segments4(k).start,...
        [ones(1,length(segments4(k).whisk_events))*cntr; ...
        ones(1,length(segments4(k).whisk_events))*cntr+1],...
        'LineWidth',3,'Color','r')
end

%% plot #2

figure
cntr = 0;

inss = ins([1 4 10 13]);
% inss = ins;

cmp = colormap('bone');
cmp = [cmp; flipud(cmp)];
spl = size(cmp,1) / 4;
cmp = [cmp(spl+1:end,:); cmp(1:spl,:)];
colormap(cmp)

S = nan(1,length(inss));
for k = inss;
    cntr = cntr + 1;
    ln = length(segments2(k).inhalation_start);
    figpos = get(gcf,'Position');
    
%     colormap('hsv')
    
    phs = resp_phase2(segments2(k).start:segments2(k).end);
    S(cntr) = subplot(length(inss),1,cntr);
    imagesc(phs')
    set(S(cntr),'XTick',[],'YTick',[])
    axis off
    
    line([segments2(k).whisk_events; segments2(k).whisk_events]-segments2(k).start,...
        [ones(1,length(segments2(k).whisk_events))*0.5; ...
        ones(1,length(segments2(k).whisk_events))*1.5],...
        'LineWidth',3,'Color','r')
end

linkaxes(S,'x')

% xlim([0 55])

%% plot #3

figure
cntr = 0;

inss = ins([1 4 10 13]);
inss = ins;
inss = (ins(2:2:14));

for k = inss;
    cntr = cntr + 1;
    ln = length(segments2(k).inhalation_start);
    figpos = get(gcf,'Position');
    X = [segments2(k).inhalation_start' segments2(k).inhalation_end' ...
        segments2(k).inhalation_end' segments2(k).inhalation_start']-segments2(k).start;
    Y = repmat([cntr cntr cntr+1 cntr+1],ln,1);
    patch(X',Y',[1 1 0.5])
    
    line([segments2(k).raw_whisk_events; segments2(k).raw_whisk_events]-segments2(k).start,...
        [ones(1,length(segments2(k).raw_whisk_events))*cntr; ...
        ones(1,length(segments2(k).raw_whisk_events))*cntr+1],...
        'LineWidth',1,'Color','r')
end


%% plot #4

figure
cntr = 0;

inss = ins([1 4 10 13]);
inss = ins;
inss = (ins(2:2:14));
inss = inss([1 5 6]);

cmp = colormap('bone');
cmp = [cmp; flipud(cmp)];
spl = size(cmp,1) / 4;
cmp = [cmp(spl+1:end,:); cmp(1:spl,:)];
colormap(cmp)

S = nan(1,length(inss));
for k = inss;
    cntr = cntr + 1;
    ln = length(segments2(k).inhalation_start);
    figpos = get(gcf,'Position');
    
%     colormap('hsv')
    
    phs = resp_phase2(segments2(k).start:segments2(k).end);
    S(cntr) = subplot(length(inss),1,cntr);
    imagesc(phs')
    set(S(cntr),'XTick',[],'YTick',[])
    axis off
    
    line([segments2(k).raw_whisk_events; segments2(k).raw_whisk_events]-segments2(k).start,...
        [ones(1,length(segments2(k).raw_whisk_events))*0.5; ...
        ones(1,length(segments2(k).raw_whisk_events))*1.5],...
        'LineWidth',1,'Color','r')
end

linkaxes(S,'x')