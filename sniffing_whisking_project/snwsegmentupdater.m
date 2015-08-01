%% load

load('C:\Balazs\_analysis\SniffWhisk\Fig_cycles\seg2to1_3.mat')
segments3_old = segments3;

%% load

load('C:\Balazs\_analysis\SniffWhisk\Fig_cycles\seg1to1_1.mat')
segments3_old = segments3;

%% load recalculated segments

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

%% make new segments3

cntr = 0;
for k = 1:length(segments3_old)
    for as = 1:length(allsegments)
        if isequal(segments3_old(k).rat,allsegments(as).rat) &&...
               isequal(segments3_old(k).session,allsegments(as).session) &&...
               (segments3_old(k).end-segments3_old(k).start)/2+segments3_old(k).start > allsegments(as).start && ...
               (segments3_old(k).end-segments3_old(k).start)/2+segments3_old(k).start < allsegments(as).end
           segments3(k) = allsegments(as);
           cntr = cntr + 1
        end
    end
end
           