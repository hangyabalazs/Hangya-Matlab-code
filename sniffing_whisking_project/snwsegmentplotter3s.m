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
    fnme = [DATAPATH 'SniffWhisk\phasehist_moreprecise2\PHASE_' rrat '_' rdate '.mat'];
    load(fnme)
    
    for k = 1:length(segments)
        segments(k).start = segments(k).start * wl / sr;
        segments(k).end = segments(k).end * wl / sr;
        
%         if segments(k).start > 6.3888e+003 && segments(k).start < 6.3890e+003...
%                 && segments(k).end > 6.3902e+003 && segments(k).end < 6.3906e+003...
%                 && segments(k).frequency > 7.0967 && segments(k).frequency < 7.0969
%             keyboard
%         end

%         if segments(k).start > 177.23 && segments(k).start < 177.25...
%                 && segments(k).end > 178.39 && segments(k).end < 178.41...
%                 && segments(k).frequency > 1.6948 && segments(k).frequency < 1.6950
%             keyboard
%         end

%         if segments(k).start > 501.7 && segments(k).start < 501.9...
%                 && segments(k).end > 502.93 && segments(k).end < 502.95...
%                 && segments(k).frequency > 7.757 && segments(k).frequency < 7.759
%             keyboard
%         end
                
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
% keyboard

%% filter

segments2 = segments(1);
for k = 1:length(segments)
    if segments(k).length > 0.5 && segments(k).group == 0.5
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
tmpp = segments2(2);
for k = 6:10
    tmp = tmpp.(fld{k});
    tmpp.(fld{k}) = tmp(1:end);
end
tmpp.start = tmpp.inhalation_start(1);
segments3 = [segments3(1:3) tmpp segments3(4:end)];


%% plot

figure
cntr = 0;

% inss = ins([1 4 10 13]);
% inss = ins;

T = 1:length(segments3);
% T = 461:483;

for k = T
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
        'LineWidth',5,'Color',[153 0 255]/255)      % blue: [0 153 255]/255, green: [0 153 0]/255
end

xlim([0 0.8])

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

%% new panel 2:1

s2inx = [20 2 18 5 fliplr([36 49 53 82 91 97])];
s2inx = [20     2    18     5  9  97    91    82    53    49    36];

segments4 = segments2(s2inx);

tt = 5;
fld = fieldnames(segments4);
tmpp = segments4(tt);
for k = 8:11
    tmp = tmpp.(fld{k});
    tmpp.(fld{k}) = tmp(2:end);
end
tmpp.start = tmpp.inhalation_start(1);
tmp = tmpp.(fld{12});
tmp(tmp<tmpp.start) = [];
tmpp.(fld{12}) = tmp;
segments4(tt) = tmpp;

tt = 4;
fld = fieldnames(segments4);
tmpp = segments4(tt);
for k = 8:11
    tmp = tmpp.(fld{k});
    tmpp.(fld{k}) = tmp(2:end);
end
tmpp.start = tmpp.inhalation_start(1);
tmp = tmpp.(fld{12});
tmp(tmp<tmpp.start) = [];
tmpp.(fld{12}) = tmp;
segments4(tt) = tmpp;

tt = 1;
fld = fieldnames(segments4);
tmpp = segments4(tt);
for k = 8:11
    tmp = tmpp.(fld{k});
    tmpp.(fld{k}) = tmp(1:end-1);
end
tmpp.start = tmpp.inhalation_start(1);
tmpp.end = tmpp.exhalation_end(end);
tmp = tmpp.(fld{12});
tmp(tmp<tmpp.start|tmp>tmpp.end) = [];
tmpp.(fld{12}) = tmp;
segments4(tt) = tmpp;
segments3 = fliplr(segments4);

%% cont.

segments5 = segments4;

s2inx2 = [2 4 6 8 10 16 21 24 27 34 36 41 43 58 65 85];
segments5 = [segments5 segments2(s2inx2)];
segments3 = fliplr(segments5);

newinx = [1 24 15 13 12 2 14 3 4 5:11 16:23 25:27];
segments5 = segments5(newinx);
segments3 = fliplr(segments5);

newinx2 = [1 3 8 2 4:7 9 10 17 18 24 21 23 25 20 22 19 26 27 11:16];
segments5 = segments5(newinx2);
segments3 = fliplr(segments5);

newinx3 = [1 2 3 7 6 8 5 4 10 13 9 16 11 12 14 15 17 18 19 20 22 21 23 25 24 26 27];
segments5 = segments5(newinx3);
segments3 = fliplr(segments5);

newinx4 = [1 3 2 4 7 6 5 8 9 11 10 13 12 14 15 16 17:21 24 22 23 25:27];
segments5 = segments5(newinx4);
segments3 = fliplr(segments5);

%% 1:2

segments4 = segments2([4 5 7 10 12 13 18 20 ...
    21 22 25 26 27 28 30 31 35 48 ...
    66 71 81 103 107 108 110 120 ...
    121 125 131 140 142 146 164 170 174 180 ...
    181 183 185 186 187 192 196 200 ...
    207 209 213:217 222 225 228 230 234 ...
    244 246 250 264 267 270 276 277 287 293 294 296:298 ...
    301 304 308 309 319 326 334 335 341 344 349 350 356 ...
    373 406 413 427 439 441 448 455 457:460 ...
    469 474 476 477]);

