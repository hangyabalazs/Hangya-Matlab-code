function nbtuningcurves3_roc(I,issave)
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
resdir = [DATAPATH 'NB' fs 'tuningcurves_roc_newdata_' align '_' group  '_' normalization fs];   % results directory
 
% PSTH
[allpsth_orig_hit allspsth_orig_hit allpsth_hit allspsth_hit outcellids_hit] = ...
    poppsth(I,resdir,align,'hit',normalization,doraster,issave); %#ok<*ASGLU>
[allpsth_orig_fa allspsth_orig_fa allpsth_fa allspsth_fa outcellids_fa] = ...
    poppsth(I,resdir,align,'fa',normalization,doraster,issave); %#ok<*ASGLU>

% Stats
mx = squeeze(nanmax(allspsth_orig_hit(:,:,1201:1300),[],3));   % max 0-100 ms of smoothed PSTH, hits
mx_fa = squeeze(nanmax(allspsth_orig_fa(:,:,1201:1300),[],3));   % max 0-100 ms of smoothed PSTH, false alarms
anova_rm(mx')   % repeated measures ANOVA
anova_rm(mx_fa')

figure
plot(mx)   % hit responses, all cells
figure
bar(1:5,[nanmean(mx_fa(:)) nanmean(mx,2)'])
hold on
errorbar(1:5,[nanmean(mx_fa(:)) nanmean(mx,2)'],[nanmean(nanse(mx_fa')) nanse(mx')],'k+')

figure
bar([1:4 6:9],[nanmean(mx_fa,2)' nanmean(mx,2)'])
hold on
errorbar([1:4 6:9],[nanmean(mx_fa,2)' nanmean(mx,2)'],[nanse(mx_fa') nanse(mx')],'k+')

C = max(mx);
C = max([mx; mx_fa]);
figure
bar([1:4 6:9],[mean(mx_fa./repmat(C,4,1),2)' mean(mx./repmat(C,4,1),2)'])
hold on
errorbar([1:4 6:9],[mean(mx_fa./repmat(C,4,1),2)' mean(mx./repmat(C,4,1),2)'],...
    [nanse(transpose(mx_fa./repmat(C,4,1))) nanse(transpose(mx./repmat(C,4,1)))],'k+')

C = nanmax([mx; mx_fa]);   % normlized to allmax for final figure
figure
bar(1:5,[nanmean(nanmean(mx_fa./repmat(C,4,1),2)) nanmean(mx./repmat(C,4,1),2)'])
hold on
errorbar(1:5,[nanmean(nanmean(mx_fa./repmat(C,4,1),2)) nanmean(mx./repmat(C,4,1),2)'],...
    [nanmean(nanse(transpose(mx_fa./repmat(C,4,1)))) nanse(transpose(mx./repmat(C,4,1)))],'k+')

% Split by depth
keyboard   % this is just command window programing from here
depth = getvalue('DVpos',outcellids_hit)';
depth_range = range(depth);
mindepth = min(depth);
inxD = depth <= mindepth + depth_range / 3;
inxM = depth > mindepth + depth_range / 3 & depth <= mindepth + 2 * depth_range / 3;
inxV = depth > mindepth + 2 * depth_range / 3;

inx = inxV;
mx = squeeze(nanmax(allspsth_orig_hit(:,inx,1201:1300),[],3));   % max 0-100 ms of smoothed PSTH, hits
mx_fa = squeeze(nanmax(allspsth_orig_fa(:,inx,1201:1300),[],3));   % max 0-100 ms of smoothed PSTH, false alarms
C = max([mx; mx_fa]);
figure
bar([1:4 6:9],[mean(mx_fa./repmat(C,4,1),2)' mean(mx./repmat(C,4,1),2)'])
hold on
errorbar([1:4 6:9],[mean(mx_fa./repmat(C,4,1),2)' mean(mx./repmat(C,4,1),2)'],...
    [nanse(transpose(mx_fa./repmat(C,4,1))) nanse(transpose(mx./repmat(C,4,1)))],'k+')

% Normalized population average
lenc = sum(inx);
allmn = squeeze(nanmean(allpsth_hit(:,inx,:),2));
allse = squeeze(nanstd(allpsth_hit(:,inx,:),[],2)) / sqrt(lenc);
allmns = squeeze(nanmean(allspsth_hit(:,inx,:),2));
allses = squeeze(nanstd(allspsth_hit(:,inx,:),[],2)) / sqrt(lenc);
wn = [-1.2 1.2];   % in seconds
dt = 0.001;
time = wn(1):dt:wn(2);   % time vector

% Plot population average
lent = size(allmn,1);
clr = summer(lent);
figure
hold on
for k = 1:lent
    errorshade(time,allmn(k,:),allse(k,:),'LineWidth',2,'LineColor',clr(k,:),...
                'ShadeColor',clr(k,:))
end
xlim([time(1) time(end)])
xlim([0 0.2])
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
function [allpsth_orig allspsth_orig allpsth allspsth outcellids] = poppsth(I,resdir,align,trialtype,normalization,doraster,issave)

% Load CellBase if indices to CELLIDLIST are passed
if isnumeric(I)
    loadcb
    I = CELLIDLIST(I);
end

% Time window
wn = [-0.1 0.5];   % in seconds
dt = 0.001;
time = wn(1):dt:wn(2);   % time vector

% Call 'ultimate_psth'
% Preallocate
problem_ids = {};   % cellIDs for which PSTH failed
problem_msg = {};   % corresponding error messages
outcellids = {};   % cellIDs corresponding to PSTH matrices
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
        [psth, spsth, spsth_se, tags spt] = ultimate_psth(cellid,'trial',cevent,wn,...
            'dt',dt,'display',true,'sigma',0.01,'parts','#StimulusDuration',...
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
            outcellids = [outcellids cellid];
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
        
        % ROC analysis
        spt_lowSI = [spt{1}; spt{2}];
        spt_highSI = [spt{3}; spt{4}];
        wns = 0.025 / dt;   % 25 ms window (in data points)
        sht = 0.005 / dt;   % shift between overlapping windows (in data points)
        nmw = floor((diff(wn)/dt-wns)/sht) + 1;   % number of overlapping windows
        [ROC P SE] = deal(nan(1,nmw));
        for w = 1:nmw
            inx1 = (w - 1) * sht + 1;
            inx2 = inx1 + wns;
            [ROC(w) P(w) SE(w)] = rocarea(nansum(spt_lowSI(:,inx1:inx2),2),...
                nansum(spt_highSI(:,inx1:inx2),2),'transform','scale','bootstrap',1000);
        end
        ROCtime = (wn(1)/dt+wns/2:sht:wn(2)/dt-wns/2) * dt;   % time vector
        ROCtime = ROCtime + wns * dt / 2;   % make the window causal
%         ROC = smooth(nan2zero(ROC),'linear',5);   % NaN when all spike counts are 0 for both distributions
        
        Hroc = figure;   % plot
%         plot(ROCtime,smooth(ROC,'linear',11),'k')
        errorshade(ROCtime(ROCtime>=-2&ROCtime<=0.6),ROC(ROCtime>=-2&ROCtime<=0.6),SE(ROCtime>=-2&ROCtime<=0.6),...
            'LineColor','k','ShadeColor','k')
        Hrocp = figure;
        plot(ROCtime,P)
        
        % Save figure
        if issave
            fnm = [resdir cellidt '_ROC_' trialtype '.fig'];
            saveas(Hroc,fnm)
            fnm = [resdir cellidt '_ROC_' trialtype '.jpg'];
            saveas(Hroc,fnm)
            fnm = [resdir cellidt '_ROCP_' trialtype '.fig'];
            saveas(Hrocp,fnm)
            fnm = [resdir cellidt '_ROCP_' trialtype '.jpg'];
            saveas(Hrocp,fnm)
            fnm = [resdir cellidt '_ROC_' trialtype '.mat'];
            save(fnm,'ROCtime','ROC','P','SE')
            close(Hroc)
            close(Hrocp)
        end
        
    catch ME
        
        % Error handling
        problem_ids{end+1} = cellid;  %#ok<AGROW> % collect cellIDs resulting in error
        problem_msg{end+1} = ME.message;  %#ok<AGROW> % store corresponding error messages
        disp(['Something went wrong for cell ' cellid '.'])
        disp(ME.message)   % display error message
    end
end

% Normalization for averaging
[lent lenc lenp] = size(allspsth_orig);
switch normalization
    case 'zscore'
        allpsth = (allpsth_orig - repmat(mean(mean(allpsth_orig,1),3),[lent,1,lenp])) ...
            ./ repmat(std(mean(allpsth_orig,1),[],3),[lent,1,lenp]);   % standardize based on mean across intensities
        allspsth = (allspsth_orig - repmat(mean(mean(allspsth_orig,1),3),[lent,1,lenp])) ...
            ./ repmat(std(mean(allspsth_orig,1),[],3),[lent,1,lenp]);  % standardize based on mean across intensities
    case 'max'
        allpsth = (allpsth_orig - repmat(min(mean(allpsth_orig,1),[],3),[lent,1,lenp])) ...
            ./ (repmat(max(mean(allpsth_orig,1),[],3),[lent,1,lenp]) - ...
            repmat(min(mean(allpsth_orig,1),[],3),[lent,1,lenp]));   % normalize based on maximum of mean across intensities
        allspsth = (allspsth_orig - repmat(min(mean(allspsth_orig,1),[],3),[lent,1,lenp])) ...
            ./ (repmat(max(mean(allspsth_orig,1),[],3),[lent,1,lenp]) - ...
            repmat(min(mean(allpsth_orig,1),[],3),[lent,1,lenp]));  % normalize based on maximum of mean across intensities
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
legend(tags)

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