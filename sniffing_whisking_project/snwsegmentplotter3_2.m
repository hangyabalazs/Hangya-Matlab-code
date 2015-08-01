%%

% load('C:\Balazs\_analysis\SniffWhisk\Fig_cycles\segments3p_v3.mat')

for k = 1:length(segments3)
    ind = 0;
    cntr = 0;
    while ind == 0
        cntr = cntr + 1;
        if ismember(segments3(k).whisk_events(2),allsegments(cntr).whisk_events) && ...
                ismember(segments3(k).whisk_events(3),allsegments(cntr).whisk_events)
            segments3(k).rat = allsegments(cntr).rat;
            segments3(k).session = allsegments(cntr).session;
            ind = 1;
%             keyboard
            
            
        end
    end
end

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
        segments(k).rat = rrat;
        segments(k).session = rdate;
    end
    
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


% kk = kk + 1;
% segments3(kk) = segments2(ins(end-16));
% for k = 6:10
%     tmp = segments3(kk).(fld{k});
%     segments3(kk).(fld{k}) = tmp(4:11);
% end
% segments3(kk).start = segments3(kk).inhalation_start(1);

kk = kk + 1;
segments3(kk) = segments2(ins(end-17));
for k = 6:10
    tmp = segments3(kk).(fld{k});
    segments3(kk).(fld{k}) = tmp(1:9);
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
segments3t = segments3;
segments3tt = segments3;

%% add more

tmpp = segments4(576);
for k = 6:10
    tmp = tmpp.(fld{k});
    tmpp.(fld{k}) = tmp(2:end);
end
tmpp.start = tmpp.inhalation_start(1);
segments3t = [segments3t(1:end) tmpp];
segments3 = segments3t;

%% add more

tmpp = segments4(573);
for k = 6:10
    tmp = tmpp.(fld{k});
    tmpp.(fld{k}) = tmp(2:end);
end
tmpp.start = tmpp.inhalation_start(1);
segments3t = [segments3t(1:end) tmpp];
segments3 = segments3t;

%% add more

tmpp = segments4(105);
for k = 6:10
    tmp = tmpp.(fld{k});
    tmpp.(fld{k}) = tmp(2:end);
end
tmpp.start = tmpp.inhalation_start(1);
segments3t = [tmpp segments3t(1:end)];
segments3 = segments3t;


%% add more

tmpp = segments4(103);
for k = 6:10
    tmp = tmpp.(fld{k});
    tmpp.(fld{k}) = tmp(1:end);
end
tmpp.start = tmpp.inhalation_start(1);
segments3t = [tmpp segments3t(1:end)];
segments3 = segments3t;

%% add more

tmpp = segments4(182);
for k = 6:10
    tmp = tmpp.(fld{k});
    tmpp.(fld{k}) = tmp(3:end);
end
tmpp.start = tmpp.inhalation_start(1);
segments3t = [segments3t(1) tmpp segments3t(2:end)];
segments3 = segments3t;

%% add more

tmpp = segments4(300);
for k = 6:10
    tmp = tmpp.(fld{k});
    tmpp.(fld{k}) = tmp(2:end);
end
tmpp.start = tmpp.inhalation_start(1);
segments3t = [segments3t(1:2) tmpp segments3t(3:end)];
segments3 = segments3t;

%% add more

tmpp = segments4(203);
for k = 6:10
    tmp = tmpp.(fld{k});
    tmpp.(fld{k}) = tmp(1:end);
end
tmpp.start = tmpp.inhalation_start(1);
segments3t = [segments3t(1:2) tmpp segments3t(3:end)];
segments3 = segments3t;

%% add more

tmpp = segments4(214);
for k = 6:10
    tmp = tmpp.(fld{k});
    tmpp.(fld{k}) = tmp(3:end);
end
tmpp.start = tmpp.inhalation_start(1);
segments3t = [segments3t(1:3) tmpp segments3t(4:end)];
segments3 = segments3t;

%% add more

tmpp = segments4(233);
for k = 6:10
    tmp = tmpp.(fld{k});
    tmpp.(fld{k}) = tmp(3:end);
end
tmpp.start = tmpp.inhalation_start(1);
segments3t = [segments3t(1:3) tmpp segments3t(4:end)];
segments3 = segments3t;

%% add more

tmpp = segments4(217);
for k = 6:10
    tmp = tmpp.(fld{k});
    tmpp.(fld{k}) = tmp(1:end);
end
tmpp.start = tmpp.inhalation_start(1);
segments3t = [segments3t(1:3) tmpp segments3t(4:end)];
segments3 = segments3t;

%% add more

tmpp = segments4(360);
for k = 6:10
    tmp = tmpp.(fld{k});
    tmpp.(fld{k}) = tmp(1:end);
end
tmpp.start = tmpp.inhalation_start(1);
segments3t = [segments3t(1:8) tmpp segments3t(9:end)];
segments3 = segments3t;

%% add more

tmpp = segments4(356);
for k = 6:10
    tmp = tmpp.(fld{k});
    tmpp.(fld{k}) = tmp(1:end);
end
tmpp.start = tmpp.inhalation_start(1);
segments3t = [segments3t(1:7) tmpp segments3t(8:end)];
segments3 = segments3t;

%% add more

tmpp = segments4(498);
for k = 6:10
    tmp = tmpp.(fld{k});
    tmpp.(fld{k}) = tmp(1:end);
end
tmpp.start = tmpp.inhalation_start(1);
segments3t = [segments3t(1:10) tmpp segments3t(11:end)];
segments3 = segments3t;

%% add more

tmpp = segments4(549);
for k = 6:10
    tmp = tmpp.(fld{k});
    tmpp.(fld{k}) = tmp(1:end);
end
tmpp.start = tmpp.inhalation_start(1);
segments3t = [segments3t(1:10) tmpp segments3t(11:end)];
segments3 = segments3t;

%% add more

tmpp = segments4(543);
for k = 6:10
    tmp = tmpp.(fld{k});
    tmpp.(fld{k}) = tmp(1:end);
end
tmpp.start = tmpp.inhalation_start(1);
segments3t = [segments3t(1:10) tmpp segments3t(11:end)];
segments3 = segments3t;

%% filter #2

segments4 = segments(1);
for k = 1:length(segments)
    if segments(k).length > 1 && segments(k).group == 1
        segments4(end+1) = segments(k);
    end
end
segments4 = segments4(2:end);


%% sort #2

lense = length(segments4);
fr = nan(1,lense);
for k = 1:lense
   fr(k) = segments4(k).frequency;
end
[sr ins] = sort(fr);
segments4 = segments4(ins);


%% plot

figure
cntr = 0;

% inss = ins([1 4 10 13]);
% inss = ins;

T = 1:length(segments3);


for k = T
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

global DATAPATH

figure
cntr = 0;
segments2 = fliplr(segments3);

% inss = ins([1 4 10 13]);
% inss = ins;
T = 1:length(segments3);

cmp = colormap('bone');
cmp = [cmp; flipud(cmp)];
spl = size(cmp,1) / 4;
cmp = [cmp(spl+1:end,:); cmp(1:spl,:)];
colormap(cmp)

S = nan(1,length(segments3));
for k = T
    cntr = cntr + 1;
    ln = length(segments2(k).inhalation_start);
    figpos = get(gcf,'Position');
    
    % Load data for phase calculation
    fnme = [DATAPATH 'SniffWhisk\phasehist8\PHASE_' segments2(k).rat '_' segments2(k).session '.mat'];
    load(fnme)
    
%     colormap('hsv')
    
    phs = resp_phase2(round(segments2(k).start*sr/wl):round(segments2(k).end*sr/wl));
    S(cntr) = subplot(length(T),1,cntr);
    imagesc(linspace(0,length(phs)/sr*wl,length(phs)),1,phs')
%     imagesc(phs')
    set(S(cntr),'XTick',[],'YTick',[])
    axis off
    wl
    
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

