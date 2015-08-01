function nbrtcurves(I,issave)
%NBRTCURVES   PSTHs for different reaction times.
%   NBRTCURVES(I,ISSAVE) calculates non-adaptive PSTHs and raster plots for
%   a set of cells aligned to 'Go' and 'No-go' tone/response onset. The
%   trials are partitioned based on reaction time quartiles. The PSTHs,
%   rasters and population average PSTHs are plotted and saved. Input
%   parameters:
%       I - index set to CELLIDLIST (see CellBase documentation); if empty
%           or not specified, a predefined set of basal forebrain cells is
%           used.
%       ISSAVE - controls saving
%
%   See also NBPERFORMANCECURVES, NBTUNINGCURVES, ULTIMATE_PSTH and
%   NBRESPONSESORTER.

%   Edit log: BH 9/19/12
 
% Pass the control to the user in case of error
dbstop if error
choosecb('NB')    % swhitch to 'NB' CellBase

% Input argument check
error(nargchk(0,2,nargin))
if nargin < 2
    issave = false;   % saving is controled by the second input argument
end
if nargin < 1
    I = [];   % initialize list of cell IDs
end
 
% List of cellIDs
group = [];
if isempty(I)
    group = 'allChAT';
    switch group
        case 'ChAT'
            selstr = ['"ChAT+"==1&"ID_PC">20&"Lr_PC"<0.15&' ...
                'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})'];
            ChAT = selectcell(selstr);   % cellIDs of ChAT+ cells
            I = ChAT;
            
        case 'pChAT'
            pChAT = {'n029_120210a_3.3' ...  % 'n029_120214b_2.7': same as one of the ChAT+ cells; 'n029_120215a_2.2': one of the ChAT+ cells recorded a day later
                'n029_120215a_3.4' 'n029_120221b_6.1' 'n029_120222b_4.1'};   % cellIDs of putative ChAT+ cells
            I = pChAT;
            
        case 'allChAT'
            selstr = ['"ChAT+"==1&"ID_PC">20&"Lr_PC"<0.15&' ...
                'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})'];
            ChAT = selectcell(selstr);   % cellIDs of ChAT+ cells
            pChAT = {'n029_120210a_3.3' ...
                'n029_120215a_3.4' 'n029_120221b_6.1' 'n029_120222b_4.1'};   % cellIDs of putative ChAT+ cells
            I = [ChAT pChAT];
            
        case 'PV'
            selstr = ['"PV+"==1&"ID_PC">20&"Lr_PC"<0.15&' ...
                'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})'];
            PV = selectcell(selstr);   % cellIDs of PV+ cells
            I = PV;
            
        case 'activated'
            activated = getgroups;   % load groups of activated and inhibited NB cells
            I = activated;
            
        case 'inhibited'
            [activated inhibited] = getgroups;   % load groups of activated and inhibited NB cells
            I = inhibited;
    end
end
I = I(:)';   % convert to row vector

% Control 'poppsth' behavior
doraster = true;
align = 'tone';  % controls alignment: 'tone' or 'response'
normalization = 'max';  % controls normalization: 'zscore' or 'max'

% Directories
global DATAPATH
fs = filesep;
resdir = [DATAPATH 'NB' fs 'rtcurves2_' align '_' group fs];   % results directory
 
% PSTH
poppsth(I,resdir,align,normalization,doraster,issave); %#ok<*ASGLU>

% -------------------------------------------------------------------------
function poppsth(I,resdir,align,normalization,doraster,issave)

% Load CellBase if indices to CELLIDLIST are passed
if isnumeric(I)
    loadcb
    I = CELLIDLIST(I);
end

% Time window
wn = [-1.2 1.2];   % in seconds
dt = 0.01;
time = wn(1):dt:wn(2);   % time vector

% Call 'ultimate_psth'
% Preallocate
problem_ids = {};   % cellIDs for which PSTH failed
problem_msg = {};   % corresponding error messages
allpsth_orig = [];  % all PSTHs, for Hit and Miss
allspsth_orig = []; % all smoothed PSTHs, for Hit and Miss
allpsth_orig2 = [];  % all PSTHs, for FA and CR
allspsth_orig2 = []; % all smoothed PSTHs, for FA and CR
NumCell = length(I);  % number of cells
for k = 1:NumCell
    cellid = I{k};
    disp(cellid)
    
    % Control PSTH and raster alignment
    switch align
        case 'tone'
            hitevent = 'StimulusOn';
            faevent = 'StimulusOn';
            sevent = 'StimulusOff';
        case 'response'
            hitevent = 'LeftWaterValveOn';
            faevent = 'LeftPortIn';
            sevent = 'StimulusOn';
    end
    
    try
        
        % Calcualte PSTH for Hits
        [psth1, spsth1, spsth_se1] = ultimate_psth(cellid,'trial',hitevent,wn,...
            'dt',dt,'display',true,'sigma',0.02,'parts','all',...
            'isadaptive',0,...
            'event_filter','selectGoRT','filterinput',[0 0.25],...
            'maxtrialno',Inf,'baselinewin',[-0.5 0],'testwin',[0 0.5],...
            'relative_threshold',0.1);
        [psth2, spsth2, spsth_se2] = ultimate_psth(cellid,'trial',hitevent,wn,...
            'dt',dt,'display',true,'sigma',0.02,'parts','all',...
            'isadaptive',0,...
            'event_filter','selectGoRT','filterinput',[0.25 0.5],...
            'maxtrialno',Inf,'baselinewin',[-0.5 0],'testwin',[0 0.5],...
            'relative_threshold',0.1);
        [psth3, spsth3, spsth_se3] = ultimate_psth(cellid,'trial',hitevent,wn,...
            'dt',dt,'display',true,'sigma',0.02,'parts','all',...
            'isadaptive',0,...
            'event_filter','selectGoRT','filterinput',[0.5 0.75],...
            'maxtrialno',Inf,'baselinewin',[-0.5 0],'testwin',[0 0.5],...
            'relative_threshold',0.1);
        [psth4, spsth4, spsth_se4] = ultimate_psth(cellid,'trial',hitevent,wn,...
            'dt',dt,'display',true,'sigma',0.02,'parts','all',...
            'isadaptive',0,...
            'event_filter','selectGoRT','filterinput',[0.75 1],...
            'maxtrialno',Inf,'baselinewin',[-0.5 0],'testwin',[0 0.5],...
            'relative_threshold',0.1);
        psth = cat(1,psth1,psth2,psth3,psth4);
        spsth = cat(1,spsth1,spsth2,spsth3,spsth4);
        spsth_se = cat(1,spsth_se1,spsth_se2,spsth_se3,spsth_se4);
        tags = {'GoRT 1st quartile','GoRT 2nd quartile','GoRT 3rd quartile','GoRT 4th quartile'};
        
        % Concatenate data from different cells
        tagsn = extractnum(tags);  % get sound intensities from tag strings
        [st stia] = sort(tagsn);
        if isempty(allpsth_orig)
            linx = 1;  % next index
        else
            linx = size(allpsth_orig,2) + 1;
        end
        for tgs = 1:length(stia)
            allpsth_orig(tgs,linx,1:length(time)) = psth(stia(tgs),:); %#ok<AGROW>
            allspsth_orig(tgs,linx,1:length(time)) = spsth(stia(tgs),:); %#ok<AGROW>
        end
        
        % Plot PSTH
        figure
        hold on
        NumPsth = size(spsth,1);
        clr = summer(NumPsth);   % colormap
        for pk = 1:NumPsth
%             plot(spsth(pk,:),'Color',clr(pk,:),'LineWidth',2)
            errorshade(time,spsth(pk,:),spsth_se(pk,:),'LineWidth',2,'LineColor',clr(pk,:),...
                'ShadeColor',clr(pk,:))
        end
        legend(tags)
%         xlim([-0.2 0.2])
        
        % Save PSTH figure
        H = gcf;
        if issave
            cellidt = regexprep(cellid,'\.','_');
            title([regexprep(cellidt,'_',' ') ' Hit'])
            fnm = [resdir cellidt '_PSTH_Hit.fig'];
            saveas(H,fnm)
        end
        close(H)
        
        % Raster plot
        if doraster
            H = figure;
            viewcell2b(cellid,'TriggerName',hitevent,'SortEvent',sevent,...
                'eventtype','behav','ShowEvents',{{sevent}},'PSTHstd','off',...
                'Partitions','#StimulusDuration&Hit','window',[-5 5])
            maximize_figure(H)
            
            % Save raster plot figure
            if issave
                fnm = [resdir cellidt '_raster_Hit.fig'];
                saveas(H,fnm)
                fnm = [resdir cellidt '_raster_Hit.jpg'];
                saveas(H,fnm)
            end
            close(H)
        end
        
    catch ME
        
        % Error handling
        problem_ids{end+1} = cellid;  %#ok<AGROW> % collect cellIDs resulting in error
        problem_msg{end+1} = ME.message;  %#ok<AGROW> % store corresponding error messages
        disp(['Something went wrong for cell ' cellid '.'])
        disp(ME.message)   % display error message
    end
end

% Normalization for averaging, Hit & Miss
% keyboard
[lent lenc lenp] = size(allspsth_orig);
switch normalization
    case 'zscore'
        allpsth = (allpsth_orig - repmat(mean(mean(allpsth_orig,1),3),[lent,1,lenp])) ...
            ./ repmat(std(mean(allpsth_orig,1),[],3),[lent,1,lenp]);   % standardize based on mean across intensities
        allspsth = (allspsth_orig - repmat(mean(mean(allspsth_orig,1),3),[lent,1,lenp])) ...
            ./ repmat(std(mean(allspsth_orig,1),[],3),[lent,1,lenp]);  % standardize based on mean across intensities
    case 'max'
        allpsth = (allpsth_orig - repmat(mean(mean(allpsth_orig,1),3),[lent,1,lenp])) ...
            ./ repmat(max(mean(allpsth_orig,1),[],3),[lent,1,lenp]);   % normalize based on maximum of mean across intensities
        allspsth = (allspsth_orig - repmat(mean(mean(allspsth_orig,1),3),[lent,1,lenp])) ...
            ./ repmat(max(mean(allspsth_orig,1),[],3),[lent,1,lenp]);  % normalize based on maximum of mean across intensities
end

% Normalized population average
allmn = squeeze(nanmean(allspsth,2));
allse = squeeze(nanstd(allspsth,[],2)) / sqrt(lenc);

% Plot population average
clr = summer(lent);
figure
hold on
for k = 1:lent
    errorshade(time,allmn(k,:),allse(k,:),'LineWidth',2,'LineColor',clr(k,:),...
                'ShadeColor',clr(k,:))
end
xlim([time(1) time(end)])
legend(tags)

keyboard

% -------------------------------------------------------------------------
function tagsn = extractnum(tags)

NumTags = length(tags);
tagsn = nan(1,NumTags);
for k = 1:NumTags
    intoken = tags{k};
    tkn = regexprep(intoken,'\s*(.*)\s*(.=)\s*(.*)','$3');   % extract stimulus intensity from the tag string
    tagsn(k) = str2double(tkn);
end

% -------------------------------------------------------------------------
function [activated inhibited] = getgroups

% Load data
global DATAPATH
load([DATAPATH 'NB\responsesorter3\allPSTH.mat'])

% PSTH statistics
activation_peak = [allstats.activation_peak];   % peak time of activation
Wpa = [allstats.Wpa];   % Mann-Whitney test for significant activation
inhibition_peak = [allstats.inhibition_peak];   % peak time of inhibition
Wpi = [allstats.Wpi];   % Mann-Whitney test for significant inhibition

% Groups of activated and inhibited cells
inx_act = Wpa < 0.01;   % significant activation
inx_inh = Wpi < 0.01;   % significant inhibition
activated = find(inx_act&~inx_inh);    % indices for activated cells
inhibited = find(inx_inh&~inx_act);    % indices for inhibited cells
ai = find(inx_act&inx_inh);
inx = activation_peak(ai) < inhibition_peak(ai);
activated_inhibited = sort(ai(inx));
inhibited_activated = ai(~inx);
activated = [activated'; activated_inhibited'];   % categorize based on first significant effect
inhibited = [inhibited'; inhibited_activated'];
activated = tags(activated);  % return cellIDs
inhibited = tags(inhibited);