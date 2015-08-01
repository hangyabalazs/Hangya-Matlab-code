function snwmodeswitch
%SNWMODESWITCH   Coupling mode statistics.
%   SNWMODESWITCH compares sniffing frequency during differnet modes of
%   sniffing-whisking phase locking (see SNWSNIFFPLOT2). It calculates
%   statistics restricted to mode switches as well.
%
%   See also SNWSNIFFPLOT2 and SNWMODESWITCH2.

% Sampling rate
sr = 1000;

% Import session data from Excel
global DATAPATH
tblfile = [DATAPATH 'SniffWhisk\data_selection.xlsx'];
[tbl0 tbl] = xlsread(tblfile);   % read session information

% Transform time variables in segments to seconds and concatenate all segments
numSessions = size(tbl,1);
allsegments = [];
sw1to2 = [];
sw2to1 = [];
sw1to05 = [];
sw05to1 = [];
sw2to05 = [];
sw05to2 = [];
gr05 = [];
gr1 = [];
gr2 = [];
for iS = 1:numSessions   % loop through all sessions
    rrat = tbl{iS,1};   % rat name
    rdate = tbl{iS,2};   % session date
    wl = tbl0(iS,3);   % constant for downsampling while segments construction - used for time unit conversion
    
    % Load segments
    fnme = [DATAPATH 'SniffWhisk\phasehist_moreprecise2\PHASE_' rrat '_' rdate '.mat'];  % load segments from the latest folder
    load(fnme)
    
    % Convert all time variables to seconds
    for k = 1:length(segments)
        segments(k).start = segments(k).start * wl / sr; %#ok<AGROW> % segment start time
        segments(k).end = segments(k).end * wl / sr; %#ok<AGROW> % segment end time
        segments(k).inhalation_start = segments(k).inhalation_start * wl / sr; %#ok<AGROW> % inhalation start times
        segments(k).inhalation_end = segments(k).inhalation_end * wl / sr; %#ok<AGROW> % inhalation end times
        segments(k).exhalation_start = segments(k).exhalation_start * wl / sr; %#ok<AGROW> % exhalation start times
        segments(k).exhalation_end = segments(k).exhalation_end * wl / sr; %#ok<AGROW> % exhalation end times
        segments(k).whisk_events = segments(k).whisk_events * wl / sr; %#ok<AGROW> % whisk event times discriminated from whisking RMS
        segments(k).raw_whisk_events = segments(k).raw_whisk_events * wl / sr; %#ok<AGROW> % whisk event times discriminated from raw data
    end
    
    % Select based on mode
    gr = [segments.group];   % group codes
    fr = [segments.frequency];   % frequencies
    gr05inx = find(gr==0.5);   % indices for 1:2 mode
    gr1inx = find(gr==1);   % indices for 1:1 mode
    gr2inx = find(gr==2);   % indices for 2:1 mode
    
    % Find mode switches
    grc = gr(1:end-1);
    grnext = gr(2:end);
    lsw1to2 = find(grc==1&grnext==2);  % 1-to-2 switches
    lsw2to1 = find(grc==2&grnext==1);  % 2-to-1 switches
    lsw1to05 = find(grc==1&grnext==0.5);  % 1-to-0.5 switches
    lsw05to1 = find(grc==0.5&grnext==1);  % 0.5-to-1 switches
    lsw2to05 = find(grc==2&grnext==0.5);  % 2-to-0.5 switches
    lsw05to2 = find(grc==0.5&grnext==2);  % 0.5-to-2 switches
    
    % Concatenate session variables
    sw1to2 = [sw1to2 length(allsegments)+lsw1to2];   %#ok<AGROW> % all 1-to-2 switches
    sw2to1 = [sw2to1 length(allsegments)+lsw2to1];   %#ok<AGROW> % all 2-to-1 switches
    sw1to05 = [sw1to05 length(allsegments)+lsw1to05];   %#ok<AGROW> % all 1-to-0.5 switches
    sw05to1 = [sw05to1 length(allsegments)+lsw05to1];   %#ok<AGROW> % all 0.5-to-1 switches
    sw2to05 = [sw2to05 length(allsegments)+lsw2to05];   %#ok<AGROW> % all 2-to-0.5 switches
    sw05to2 = [sw05to2 length(allsegments)+lsw05to2];   %#ok<AGROW> % all 0.5-to-2 switches
    gr05 = [gr05 length(allsegments)+gr05inx];   %#ok<AGROW>  % all indices for 1:2 mode
    gr1 = [gr1 length(allsegments)+gr1inx];   %#ok<AGROW>  % all indices for 1:1 mode
    gr2 = [gr2 length(allsegments)+gr2inx];   %#ok<AGROW>  % all indices for 2:1 mode
    allsegments = [allsegments segments];   %#ok<AGROW> % all segments from all data sessions
end

% Replace segments variable with all recorded segments
segments = allsegments;
gr = [segments.group];   % group codes corresponding to all segments
fr = [segments.frequency];   % frequencies corresponding to all segments

% Mode frequency statistics
figure
boxplot([fr(gr05) fr(gr1) fr(gr2)],[zeros(size(gr05)) ones(size(gr1)) ones(size(gr2))*2],'labels',[{'1:2'} {'1:1'} {'2:1'}]);

% Plot with SD
figure
hold on
B05 = bar((1:2),[mean(fr(gr05)) mean(fr(gr05)/2)]);
B1 = bar((4:5),[mean(fr(gr1)) mean(fr(gr1))]);
B2 = bar((7:8),[mean(fr(gr2)) mean(fr(gr2)*2)]);
set(B05,'FaceColor','w','EdgeColor',[153 0 255]/255,'LineWidth',2)
set(B1,'FaceColor','w','EdgeColor',[0 153 0]/255,'LineWidth',2)
set(B2,'FaceColor','w','EdgeColor',[0 153 255]/255,'LineWidth',2)

E05 = errorbar((1:2),[mean(fr(gr05)) mean(fr(gr05)/2)],...
    [std(fr(gr05)) std(fr(gr05)/2)],...
    '+','Color',[153 0 255]/255,'LineWidth',2);
errorbar_tick(E05,0)
e05 = get(E05,'Children');
set(e05(2),'LineStyle',':')
E1 = errorbar((4:5),[mean(fr(gr1)) mean(fr(gr1))],...
    [std(fr(gr1)) std(fr(gr1))],...
    '+','Color',[0 153 0]/255,'LineWidth',2);
errorbar_tick(E1,0)
e1 = get(E1,'Children');
set(e1(2),'LineStyle',':')
E2 = errorbar((7:8),[mean(fr(gr2)) mean(fr(gr2)*2)],...
    [std(fr(gr2)) std(fr(gr2)*2)],...
    '+','Color',[0 153 255]/255,'LineWidth',2);
errorbar_tick(E2,0)
e2 = get(E2,'Children');
set(e2(2),'LineStyle',':')

% Plot with SE
figure
hold on
B05 = bar((1:2),[mean(fr(gr05)) mean(fr(gr05)/2)]);
B1 = bar((4:5),[mean(fr(gr1)) mean(fr(gr1))]);
B2 = bar((7:8),[mean(fr(gr2)) mean(fr(gr2)*2)]);
set(B05,'FaceColor','w','EdgeColor',[153 0 255]/255,'LineWidth',2)
set(B1,'FaceColor','w','EdgeColor',[0 153 0]/255,'LineWidth',2)
set(B2,'FaceColor','w','EdgeColor',[0 153 255]/255,'LineWidth',2)

E05 = errorbar((1:2),[mean(fr(gr05)) mean(fr(gr05)/2)],...
    [std(fr(gr05))/sqrt(length(gr05)) std(fr(gr05)/2)/sqrt(length(gr05))],...
    '+','Color',[153 0 255]/255,'LineWidth',2);
errorbar_tick(E05,0)
E1 = errorbar((4:5),[mean(fr(gr1)) mean(fr(gr1))],...
    [std(fr(gr1))/sqrt(length(gr1)) std(fr(gr1))/sqrt(length(gr1))],...
    '+','Color',[0 153 0]/255,'LineWidth',2);
errorbar_tick(E1,0)
E2 = errorbar((7:8),[mean(fr(gr2)) mean(fr(gr2)*2)],...
    [std(fr(gr2))/sqrt(length(gr2)) std(fr(gr2)*2)/sqrt(length(gr2))],...
    '+','Color',[0 153 255]/255,'LineWidth',2);
errorbar_tick(E2,0)

% Comparison of sniffing frequency before and after the switches
boxstat(fr(sw1to2),fr(sw1to2+1),'before switch','after switch')

% Sniffing frequency distribution of the modes
edges = 2:0.5:15;
cnts = (edges(1:end-1) + edges(2:end)) / 2;
b05 = histc(fr(gr05),edges);  % frequency histogram for 1:2
b05 = b05(1:end-1);
b1 = histc(fr(gr1),edges+0.15);  % frequency histogram for 1:1
b1 = b1(1:end-1);
b2 = histc(fr(gr2),edges+0.3);  % frequency histogram for 2:1
b2 = b2(1:end-1);
figure
hold on
stairs(cnts,b05/sum(b05),'Color',[153 0 255]/255,'LineWidth',2)  % plot
stairs(cnts+0.15,b1/sum(b1),'Color',[0 153 0]/255,'LineWidth',2)
stairs(cnts+0.3,b2/sum(b2),'Color',[0 153 255]/255,'LineWidth',2)