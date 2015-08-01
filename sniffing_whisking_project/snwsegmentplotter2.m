%% load

% Import session data from Excel
global DATAPATH
tblfile = [DATAPATH 'SniffWhisk\data_selection.xlsx'];
[tbl0 tbl] = xlsread(tblfile);

% Call 'snwdisc'
numSessions = size(tbl,1);
allsegments = [];
for iS = 1:numSessions
    rrat = tbl{iS,1};
    rdate = tbl{iS,2};
    
    % Load data for phase calculation
    fnme = [DATAPATH 'SniffWhisk\phasehist3\PHASE_' rrat '_' rdate '.mat'];
    load(fnme)
    
    allsegments = [allsegments segments];
end
segments = allsegments;

%% filter

segments2 = segments(1);
for k = 1:length(segments)
    if segments(k).length > 2 && segments(k).group == 1
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

%% restrict data

segments3 = segments2(ins(end));
fld = fieldnames(segments3);
for k = 6:10
    tmp = segments3(1).(fld{k});
    segments3(1).(fld{k}) = tmp(14:22);
end
segments3(1).start = segments3(1).inhalation_start(1);

segments3(2) = segments2(ins(end-1));
for k = 6:10
    tmp = segments3(2).(fld{k});
    segments3(2).(fld{k}) = tmp(10:22);
end
segments3(2).start = segments3(2).inhalation_start(1);

tmpp = segments3(1);
segments3(1) = segments3(2);
segments3(2) = tmpp;

kk = 3;
segments3(kk) = segments2(ins(end-kk+1));
for k = 6:10
    tmp = segments3(kk).(fld{k});
    segments3(kk).(fld{k}) = tmp(19:26);
end
segments3(kk).start = segments3(kk).inhalation_start(1);

kk = 4;
segments3(kk) = segments2(ins(end-kk+1));
for k = 6:10
    tmp = segments3(kk).(fld{k});
    segments3(kk).(fld{k}) = tmp(3:11);
end
segments3(kk).start = segments3(kk).inhalation_start(1);

%% restrict data #2

kk = 0;

kk = kk + 1;
segments3(kk) = segments2(ins(end-16));
for k = 6:10
    tmp = segments3(kk).(fld{k});
    segments3(kk).(fld{k}) = tmp(4:11);
end
segments3(kk).start = segments3(kk).inhalation_start(1);

kk = kk + 1;
segments3(kk) = segments2(ins(end-17));
for k = 6:10
    tmp = segments3(kk).(fld{k});
    segments3(kk).(fld{k}) = tmp(1:9);
end
segments3(kk).start = segments3(kk).inhalation_start(1);

% kk = kk + 1;
% segments3(kk) = segments2(ins(end));
% for k = 6:10
%     tmp = segments3(kk).(fld{k});
%     segments3(kk).(fld{k}) = tmp(12:23);
% end
% segments3(kk).start = segments3(kk).inhalation_start(1);

kk = kk + 1;
segments3(kk) = segments2(ins(end-1));
for k = 6:10
    tmp = segments3(kk).(fld{k});
    segments3(kk).(fld{k}) = tmp(11:22);
end
segments3(kk).start = segments3(kk).inhalation_start(1);

kk = kk + 1;
segments3(kk) = segments2(ins(end-2));
for k = 6:10
    tmp = segments3(kk).(fld{k});
    segments3(kk).(fld{k}) = tmp(16:24);
end
segments3(kk).start = segments3(kk).inhalation_start(1);

kk = kk + 1;
segments3(kk) = segments2(ins(end-5));
for k = 6:10
    tmp = segments3(kk).(fld{k});
    segments3(kk).(fld{k}) = tmp(5:10);
end
segments3(kk).start = segments3(kk).inhalation_start(1);

kk = kk + 1;
segments3(kk) = segments2(ins(end-6));
fld = fieldnames(segments3);
for k = 6:10
    tmp = segments3(kk).(fld{k});
    segments3(kk).(fld{k}) = tmp(18:30);
end
segments3(kk).start = segments3(kk).inhalation_start(1);

kk = kk + 1;
segments3(kk) = segments2(ins(end-15));
for k = 6:10
    tmp = segments3(kk).(fld{k});
    segments3(kk).(fld{k}) = tmp(1:5);
end
segments3(kk).start = segments3(kk).inhalation_start(1);

kk = kk + 1;
segments3(kk) = segments2(ins(end-10));
for k = 6:10
    tmp = segments3(kk).(fld{k});
    segments3(kk).(fld{k}) = tmp(5:11);
end
segments3(kk).start = segments3(kk).inhalation_start(1);

kk = kk + 1;
segments3(kk) = segments2(ins(end-13));
for k = 6:10
    tmp = segments3(kk).(fld{k});
    segments3(kk).(fld{k}) = tmp(7:16);
end
segments3(kk).start = segments3(kk).inhalation_start(1);

kk = kk + 1;
segments3(kk) = segments2(ins(end-11));
for k = 6:10
    tmp = segments3(kk).(fld{k});
    segments3(kk).(fld{k}) = tmp(28:32);
end
segments3(kk).start = segments3(kk).inhalation_start(1);


kk = kk + 1;
segments3(kk) = segments2(ins(end-18));
for k = 6:10
    tmp = segments3(kk).(fld{k});
    segments3(kk).(fld{k}) = tmp(16:20);
end
segments3(kk).start = segments3(kk).inhalation_start(1);




segments3 = fliplr(segments3);

%% plot

figure
cntr = 0;

% inss = ins([1 4 10 13]);
% inss = ins;

for k = 1:length(segments3);
    cntr = cntr + 1;
    ln = length(segments3(k).inhalation_start);
    figpos = get(gcf,'Position');
    X = [segments3(k).inhalation_start' segments3(k).inhalation_end' ...
        segments3(k).inhalation_end' segments3(k).inhalation_start']-segments3(k).start;
    Y = repmat([cntr cntr cntr+1 cntr+1],ln,1);
    patch(X',Y',[1 1 0.5])
    
    line([segments3(k).whisk_events; segments3(k).whisk_events]-segments3(k).start,...
        [ones(1,length(segments3(k).whisk_events))*cntr; ...
        ones(1,length(segments3(k).whisk_events))*cntr+1],...
        'LineWidth',3,'Color','r')
end

% xlim([0 55])

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