%% load

load('C:\Balazs\_analysis\SniffWhisk\phasehist3\PHASE_R3_2008-07-09_14-46-47_R3.mat')

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

%% plot

figure
cntr = 0;

% inss = ins([1 4 10 13]);
inss = ins;

for k = inss;
    cntr = cntr + 1;
    ln = length(segments2(k).inhalation_start);
    figpos = get(gcf,'Position');
    X = [segments2(k).inhalation_start' segments2(k).inhalation_end' ...
        segments2(k).inhalation_end' segments2(k).inhalation_start']-segments2(k).start;
    Y = repmat([cntr cntr cntr+1 cntr+1],ln,1);
    patch(X',Y',[1 1 0.5])
    
    line([segments2(k).whisk_events; segments2(k).whisk_events]-segments2(k).start,...
        [ones(1,length(segments2(k).whisk_events))*cntr; ...
        ones(1,length(segments2(k).whisk_events))*cntr+1],...
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