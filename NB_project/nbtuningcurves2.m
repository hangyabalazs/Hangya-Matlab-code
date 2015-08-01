function nbtuningcurves2(I,issave)
%NBTUNINGCURVES   PSTHs for different tone intensities.
%   NBTUNINGCURVES(I,ISSAVE) calculates non-adaptive PSTHs and raster plots
%   for a set of cells aligned to 'Go' and 'No-go' tone/response onset
%   ('Hit' and 'False Alarm'). The trials are partitioned based on stimulus
%   intensity. The PSTHs, rasters and population average PSTHs are plotted
%   and saved.
%   Input parameters: 
%       I - index set to CELLIDLIST (see CellBase documentation); if empty
%           or not specified, a predefined set of basal forebrain cells is
%           used.
%       ISSAVE - controls saving
%
%   See also ULTIMATE_PSTH and NBRESPONSESORTER.

%   Edit log: BH 9/5/12, 2/13/14
 
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
            selstr = ['"ChAT+"==1&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
                'ismember("session_type",{''gonogo'',''feedbackdelay''})&' ...
                'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})'];
            ChAT = selectcell(selstr);   % cellIDs of ChAT+ cells
            ChAT = [ChAT 'n045_121217x_4.6'];   % clustered based on light-evoked spikes
            I = ChAT;
            
        case 'pChAT'
            selstr = ['"pChAT+"==1&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
                'ismember("session_type",{''gonogo'',''feedbackdelay''})&' ...
                'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})'];
            pChAT = selectcell(selstr);  % putative
            I = pChAT;
            
        case 'allChAT'
            selstr = ['"ChAT+"==1&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
                'ismember("session_type",{''gonogo'',''feedbackdelay''})&' ...
                'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})'];
            ChAT = selectcell(selstr);   % cellIDs of ChAT+ cells
            ChAT = [ChAT 'n045_121217x_4.6'];   % clustered based on light-evoked spikes
            selstr = ['"pChAT+"==1&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
                'ismember("session_type",{''gonogo'',''feedbackdelay''})&' ...
                'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})'];
            pChAT = selectcell(selstr);  % putative
            I = [ChAT pChAT];
            
        case 'PV'
            selstr = ['"PV+"==1&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
                 'ismember("session_type",{''gonogo'',''feedbackdelay''})&'...
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
doraster = false;
align = 'feedback';  % controls alignment: 'tone' or 'response' or 'feedback'
normalization = 'max';  % controls normalization: 'zscore' or 'max'

% Directories
global DATAPATH
fs = filesep;
resdir = [DATAPATH 'NB' fs 'tuningcurves_new3_' align '_' group  '_' normalization fs];   % results directory
 
% PSTH
[allpsth_orig_hit allspsth_orig_hit] = poppsth(I,resdir,align,'hit',normalization,doraster,issave); %#ok<*ASGLU>
[allpsth_orig_fa allspsth_orig_fa] = poppsth(I,resdir,align,'fa',normalization,doraster,issave); %#ok<*ASGLU>
[allpsth_hit allspsth_hit allpsth_fa allspsth_fa] = ...
    normalize(allpsth_orig_hit,allpsth_orig_fa,allspsth_orig_hit,allspsth_orig_fa,...
    normalization,issave);

% Stats
mx = squeeze(max(allspsth_orig_hit(:,:,1201:1300),[],3));   % max 0-100 ms of smoothed PSTH, hits
mx_fa = squeeze(max(allspsth_orig_fa(:,:,1201:1300),[],3));   % max 0-100 ms of smoothed PSTH, false alarms
anova_rm(mx')   % repeated measures ANOVA
anova_rm(mx_fa')

figure
plot(mx)   % hit responses, all cells
figure
bar(1:5,[mean(mx_fa(:)) mean(mx,2)'])
hold on
errorbar(1:5,[mean(mx_fa(:)) mean(mx,2)'],[mean(nanse(mx_fa')) nanse(mx')],'k+')

% -------------------------------------------------------------------------
function [allpsth_orig allspsth_orig] = poppsth(I,resdir,align,trialtype,normalization,doraster,issave)

% Load CellBase if indices to CELLIDLIST are passed
if isnumeric(I)
    loadcb
    I = CELLIDLIST(I);
end

% Time window
wn = [-1.2 1.2];   % in seconds
dt = 0.001;
time = wn(1):dt:wn(2);   % time vector

% Call 'ultimate_psth'
% Preallocate
problem_ids = {};   % cellIDs for which PSTH failed
problem_msg = {};   % corresponding error messages
allpsth_orig = [];  % all PSTHs, sorted according to stimulus intensity
allspsth_orig = []; % all smoothed PSTHs, sorted according to stimulus intensity
NumCell = length(I);  % number of cells
for k = 1:NumCell
    cellid = I{k};
    disp(cellid)
    
    % Control PSTH and raster alignment
    switch align
        case 'tone'
            hitevent = 'StimulusOn'; %#ok<NASGU>
            faevent = 'StimulusOn'; %#ok<NASGU>
            sevent = 'StimulusOff';
        case 'response'
            hitevent = 'LeftWaterValveOn'; %#ok<NASGU>
            faevent = 'LeftPortIn'; %#ok<NASGU>
            sevent = 'StimulusOn';
        case 'feedback'
            hitevent = findAlignEvent_posfeedback_gonogo(cellid); %#ok<NASGU>
            faevent = findAlignEvent_negfeedback_gonogo(cellid); %#ok<NASGU>
            sevent = 'StimulusOn';
    end
    cevent = eval([trialtype 'event']);  % hitevent or faevent
    hitfilter = 'Hit==1'; %#ok<NASGU>
    fafilter = 'FalseAlarm==1'; %#ok<NASGU>
    cfilter = eval([trialtype 'filter']);   % hitfilter or fafilter
        
    try
        
        % Calcualte PSTH for Hits
        [psth, spsth, spsth_se, tags] = ultimate_psth(cellid,'trial',cevent,wn,...
            'dt',dt,'display',true,'sigma',0.02,'parts','#StimulusDuration',...
            'isadaptive',2,'event_filter','custom','filterinput',cfilter,...
            'maxtrialno',Inf,'baselinewin',[-0.5 0],'testwin',[0 0.5],...
            'relative_threshold',0.1);
        
        % Concatenate data from different cells
        if isequal(length(tags),4)
            tagsn = extractnum(tags);  % get sound intensities from tag strings
            [st stia] = sort(tagsn);
            if isempty(allpsth_orig)
                linx = 1;  % next index
            else
                linx = size(allpsth_orig,2) + 1;
            end
            for tgs = 1:4
                allpsth_orig(tgs,linx,1:length(time)) = psth(stia(tgs),:); %#ok<AGROW>
                allspsth_orig(tgs,linx,1:length(time)) = spsth(stia(tgs),:); %#ok<AGROW>
            end
        end
        
        % Plot PSTH
        figure
        hold on
        NumPsth = size(spsth,1);
        clr = summer(NumPsth);   % colormap
        for pk = 1:NumPsth
            errorshade(time,spsth(pk,:),spsth_se(pk,:),'LineWidth',2,'LineColor',clr(pk,:),...
                'ShadeColor',clr(pk,:))
        end
        legend(tags)
        
        % Save PSTH figure
        H = gcf;
        if issave
            cellidt = regexprep(cellid,'\.','_');
            title([regexprep(cellidt,'_',' ') ' ' trialtype])
            fnm = [resdir cellidt '_PSTH_' trialtype '.fig'];
            saveas(H,fnm)
        end
        close(H)
        
        % Raster plot
        if doraster
            H = figure;
            viewcell2b(cellid,'TriggerName',cevent,'SortEvent',sevent,...
                'eventtype','behav','ShowEvents',{{sevent}},'PSTHstd','off',...
                'Partitions','#StimulusDuration&Hit','window',[-5 5])
            maximize_figure(H)
            
            % Save raster plot figure
             if issave
                fnm = [resdir cellidt '_raster_' trialtype '.fig'];
                saveas(H,fnm)
                fnm = [resdir cellidt '_raster_' trialtype '.jpg'];
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

% -------------------------------------------------------------------------
function [allpsth_hit allspsth_hit allpsth_fa allspsth_fa] = ...
    normalize(allpsth_orig_hit,allpsth_orig_fa,allspsth_orig_hit,allspsth_orig_fa,...
    normalization,issave)

% Time window
wn = [-1.2 1.2];   % in seconds
dt = 0.001;
time = wn(1):dt:wn(2);   % time vector

% Normalization for averaging
allpsth_orig = [allpsth_orig_hit; allpsth_orig_fa];
allspsth_orig = [allspsth_orig_hit; allspsth_orig_fa];
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
allmn = squeeze(nanmean(allpsth,2));
allse = squeeze(nanstd(allpsth,[],2)) / sqrt(lenc);
allmns = squeeze(nanmean(allspsth,2));
allses = squeeze(nanstd(allspsth,[],2)) / sqrt(lenc);

% Plot population average
clr = summer(lent);
figure
hold on
for k = 1:lent
    errorshade(time,allmn(k,:),allse(k,:),'LineWidth',2,'LineColor',clr(k,:),...
                'ShadeColor',clr(k,:))
end
xlim([time(1) time(end)])
% legend(tags)

% figure
% hold on
% for k = 1:lent
%     plot(time,smooth(allmn(k,:),'linear',11),'LineWidth',2,'Color',clr(k,:))
% end
% xlim([time(1) time(end)])
% legend(tags)

% Save PSTH figure
H = gcf;
if issave
    fnm = [resdir 'meanPSTH_' trialtype '.fig'];
    saveas(H,fnm)
end
close(H)

% Plot population average
clr = summer(lent);
figure
hold on
for k = 1:lent
    errorshade(time,allmns(k,:),allses(k,:),'LineWidth',2,'LineColor',clr(k,:),...
        'ShadeColor',clr(k,:))
end
xlim([time(1) time(end)])
legend(tags)

% Save PSTH figure
H = gcf;
if issave
    fnm = [resdir 'meanSPSTH_' trialtype '.fig'];
    saveas(H,fnm)
end
close(H)

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