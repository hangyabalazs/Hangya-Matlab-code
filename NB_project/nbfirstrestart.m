function nbfirstrestart
%NBFIRSTRESTART   Raster plot and PSTH aligned to the first ITI restart.
%   NBFIRSTRESTART plots PSTH and raster plot aligned to the first
%   foreperiod lick if the foreperiod was restarted.
%
%   See also NBATTENTIONCELLS3 and ULTIMATE_PSTH.

%   Balazs Hangya, Cold Spring Harbor Laboratory
%   1 Bungtown Road, Cold Spring Harbor
%   balazs.cshl@gmail.com
%   24-Jul-2014

%   Edit log: BH 7/24/14

% Directories
global DATAPATH
resdir = fullfile(DATAPATH,'NB','firstrestart',filesep);   % results directory
issave = true;

% Cell types
ChAT = selectcell(['"ChAT+"==1&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})']);  % identified ChAT+ cells
ChAT = [ChAT 'n045_121217x_4.6'];   % clustered based on light-evoked spikes

pChAT = selectcell(['"pChAT+"==1&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})']);  % putative ChAT+ cells

PV = selectcell(['"PV+"==1&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})']);  % identified PV+ cells

cellids = [ChAT pChAT];   % identified and putative cholinergic neurons

% PSTH, raster plot, ROC analysis
main(cellids,resdir,issave)

% -------------------------------------------------------------------------
function [error_list1 error_list2] = main(cellids,resdir,issave)
% Based on NBSTIMRASTER2 and NBITIFR

% PSTH aligned to first restart
[allpsth_orig allpsth allspsth_orig allspsth allstats tags wn...
    error_list1] = poppsth(cellids,'FirstRestart',[],[],resdir,issave); %#ok<*ASGLU>
time = linspace(wn(1)*1000,wn(2)*1000,size(allpsth,2)); %#ok<NASGU>
fnm = [resdir 'allPSTH_firstrestart.mat'];
if issave
    save(fnm,'allpsth_orig','allpsth','allspsth_orig','allspsth','allstats','time',...
        'tags','error_list1')
end

% % PSTH aligned to stimulus on
% [allpsth_orig allpsth allspsth_orig allspsth allstats tags wn...
%     error_list2] = poppsth(cellids,'StimulusOn','LastITIBegins',[],resdir,issave); %#ok<*ASGLU>
% time = linspace(wn(1)*1000,wn(2)*1000,size(allpsth,2)); %#ok<NASGU>
% fnm = [resdir 'allPSTH_StimulusOn.mat'];
% if issave
%     save(fnm,'allpsth_orig','allpsth','allspsth_orig','allspsth','allstats','time',...
%         'tags','error_list2')
% end

% Create time-stamped error list in base workspace
dsr = datestr(now);  % date stamp
dsr = regexprep(dsr,':','_');
list_name1 = ['error_list1_' dsr];  % variable name
list_name1 = regexprep(list_name1,'-','_');
list_name1 = regexprep(list_name1,' ','_');
assignin('base',list_name1,error_list1)   % assign error list in base workspace
error_fname1 = fullfile(resdir,[list_name1 '.mat']);   % file name
save(error_fname1,'error_list1')   % save error list
% list_name2 = ['error_list2_' dsr];  % variable name
% list_name2 = regexprep(list_name2,'-','_');
% list_name2 = regexprep(list_name2,' ','_');
% assignin('base',list_name2,error_list2)   % assign error list in base workspace
% error_fname2 = fullfile(resdir,[list_name2 '.mat']);   % file name
% save(error_fname2,'error_list2')   % save error list

% -------------------------------------------------------------------------
function [allpsth allpsth2 allspsth allspsth2 allstats ...
    ok_ids wn error_list] = poppsth(cellids,alignevent,...
    first_event,last_event,resdir,issave)

% Time window
wn = [-3 3];   % in seconds
dt = 0.001;   % resolution, in seconds
time = wn(1):dt:wn(2);   % time vector
clr = [0.5 0 0; 1 0 0; 1 0.5 0; 1 1 0];

% Call 'main'
% Preallocate
[allpsth allspsth] = deal([]);
allstats = struct([]);
error_list = struct('cellid',{},'message',{},'exception',{});  % keep a list of caught errors 
errinx = 0;
ok_ids = {};
NumCells = length(cellids);
for k = 1:NumCells
    cellid = cellids{k};   % cell ID
    disp(cellid)
    
    try
        
        % Add new event for FirstRestart
        TS = loadcb(cellid,'EVENTSPIKES');   % load prealigned spikes
        if isequal(findcellstr(TS.events(:,1),'FirstRestart'),0)   % prealign spikes to trial start if it has not happened before
            prealignSpikes(cellid,'FUNdefineEventsEpochs',...
                @defineEventsEpochs_firstrestart,'filetype','behav',...
                'ifsave',1,'ifappend',1)
        end
        
        % Calcualte PSTH
        [psth, spsth, spsth_se, ~, spt, stats] = ...
            ultimate_psth(cellid,'trial',alignevent,wn,...
            'dt',dt,'display',false,'sigma',0.02,'parts','all','isadaptive',2,...
            'maxtrialno',Inf,'first_event',first_event,'last_event',last_event,...
            'baselinewin',[-0.5 0],'testwin',[0 0.5],'relative_threshold',0.1);
        figure
        errorshade(time,spsth,spsth_se,'LineWidth',2,...
            'LineColor',clr(1,:),'ShadeColor',clr(1,:))
        hold on
        set(gca,'Color',[0.8 0.8 0.8])
        
        % Concatenate data from different cells
        allpsth = [allpsth; psth]; %#ok<AGROW>
        allspsth = [allspsth; spsth]; %#ok<AGROW>
        allstats = [allstats; stats]; %#ok<AGROW>
        ok_ids{end+1} = cellid; %#ok<AGROW>
        
        % Save figure
        H = gcf;
        cellidt = regexprep(cellid,'\.','_');
        if issave
            fnm = [resdir cellidt '_' alignevent '_PSTH.fig'];
            saveas(H,fnm)
        end
        close(H)
        
        % Raster plot with sorted trials
        Hr2 = figure;
        viewcell2b(cellids(k),'TriggerName',alignevent,'SortEvent','StimulusOff',...
            'eventtype','behav','ShowEvents',{{'StimulusOn','LeftPortIn'}},...
            'Partitions','all','window',wn)
        maximize_figure(Hr2)
        cla   % replace PSTH
        errorshade(time,spsth,spsth_se,'LineWidth',2,...
            'LineColor',clr(1,:),'ShadeColor',clr(1,:))
        hold on
        set(gca,'Color',[0.8 0.8 0.8])
        axis tight
        set(gcf,'Renderer','OpenGL')
        A = findobj(allchild(gcf),'Type','axes');  % overlay zero line
        axes(A(5));
        line([0 0],ylim,'LineWidth',3,'Color',[0 0 0.75])
        fnm = [resdir cellidt '_' alignevent '_rasterPSTH.fig'];   % save fig
        saveas(Hr2,fnm)
        fnm = [resdir cellidt '_' alignevent '_rasterPSTH.jpg'];   % save jpg
        if issave
            set(Hr2,'PaperPositionMode','auto')
            set(Hr2,'InvertHardcopy','off')
            print(Hr2,'-djpeg',fnm)
        end
        close(Hr2)
        
    catch ME
        
        % Error handling
        errinx = errinx + 1;  % error counter
        error_list(errinx).cellid = cellid;   % record problematic cell ID
        error_list(errinx).message = ME.message;   % error message
        error_list(errinx).exception = ME;   % exception structure
        disp(['Something went wrong for cell ' cellid '.'])
        disp(ME.message)
    end
end

% Standardize
allpsth2 = allpsth;
allspsth2 = allspsth;
for k = 1:size(allpsth,1)
    allpsth2(k,:) = standardize(allpsth(k,:));
    allspsth2(k,:) = standardize(allspsth(k,:));
end