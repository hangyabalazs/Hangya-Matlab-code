%% tagging

% tagging
% use insertdata to add H, D_KL and R to CellBase
% check clusters and use insertdata to add validity to CellBase
% go through the rasters and identify candidates (cells.xlsx)

% taggedprop
% use insertdata to add ChAT+ (1 for primary and 2 for alternative
% sessions) and alternative_sessions to CellBase

%% raster plots

% fig_raster_caller

%% FA PSTH

addanalysis(@ultimate_psth,...
  'mandatory',{'trial' 'DeliverFeedback' [-0.6 0.6]},...
  'property_names',{'FA_psth' 'FA_psth_stats'},'output_subset',[1 6],...
  'arglist',{'dt',0.001;'display',false;'sigma',0.02;'parts','all';'isadaptive',2;...
  'event_filter','custom';'filterinput','FalseAlarm==1';'maxtrialno',Inf;...
  'baselinewin',[-0.5 0];'testwin',[0 0.5];'relative_threshold',0.1});

R = runanalysis(@ultimate_psth,...
'trial', 'DeliverFeedback', [-0.6 0.6],...
'dt',0.001,'display',false,'sigma',0.02,'parts','all','isadaptive',2,...
'event_filter','custom','filterinput','FalseAlarm==1','maxtrialno',Inf,...
'baselinewin',[-0.5 0],'testwin',[0 0.5],'relative_threshold',0.1,'cellids',CELLIDLIST(3851));

%% FA_psth, FA_psth_stats

% N = 3456:3852;
N = 338:487;
R = cell(length(N),3);
cntr = 0;
for iC = N
    cntr = cntr + 1;
    cellid = CELLIDLIST{iC};
    disp(cellid)
    if ismember(getvalue('session_type',cellid),{'feedbackdelay' 'gonogo'})
        R{cntr,1} = cellid;
        R(cntr,2:3) = runanalysis(@ultimate_psth,...
            'trial', 'DeliverFeedback', [-0.6 0.6],...
            'dt',0.001,'display',false,'sigma',0.02,'parts','all','isadaptive',2,...
            'event_filter','custom','filterinput','FalseAlarm==1','maxtrialno',Inf,...
            'baselinewin',[-0.5 0],'testwin',[0 0.5],'relative_threshold',0.1,...
            'cellids',CELLIDLIST(iC),'outputargs',[1 6]);
    end
end
vinx = cellfun(@(s)~isempty(s),R(:,1));
insertdata(R(vinx,:),'type','prop','name',{'FA_psth' 'FA_psth_stats'})

%% Hit_psth, Hit_psth_stats

% N = 3456:3852;
N = 338:487;
R = cell(length(N),3);
cntr = 0;
for iC = N
    cntr = cntr + 1;
    cellid = CELLIDLIST{iC};
    disp(cellid)
    if ismember(getvalue('session_type',cellid),{'feedbackdelay' 'gonogo'})
        R{cntr,1} = cellid;
        R(cntr,2:3) = runanalysis(@ultimate_psth,...
            'trial', 'DeliverFeedback', [-0.6 0.6],...
            'dt',0.001,'display',false,'sigma',0.02,'parts','all','isadaptive',2,...
            'event_filter','custom','filterinput','Hit==1','maxtrialno',Inf,...
            'baselinewin',[-0.5 0],'testwin',[0 0.5],'relative_threshold',0.1,...
            'cellids',CELLIDLIST(iC),'outputargs',[1 6]);
    end
end
vinx = cellfun(@(s)~isempty(s),R(:,1));
insertdata(R(vinx,:),'type','prop','name',{'Hit_psth' 'Hit_psth_stats'})

%% Hit PSTH

addanalysis(@ultimate_psth,...
  'mandatory',{'trial' 'DeliverFeedback' [-0.6 0.6]},...
  'property_names',{'Hit_psth' 'Hit_psth_stats'},'output_subset',[1 6],...
  'arglist',{'dt',0.001;'display',false;'sigma',0.02;'parts','all';'isadaptive',2;...
  'event_filter','custom';'filterinput','Hit==1';'maxtrialno',Inf;...
  'baselinewin',[-0.5 0];'testwin',[0 0.5];'relative_threshold',0.1});

%% surprise

% nbtuningcurves3