function nbsalience2(I,issave)
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
align1 = 'tone';  % controls alignment: 'tone' or 'response' or 'feedback
align2 = 'feedback';

% Directories
global DATAPATH
fs = filesep;
resdir = [DATAPATH 'NB' fs 'salience_newdata_' group fs];   % results directory
 
% PSTH
[baseline_tone maxvalue_tone] = poppsth(I,resdir,align1,'hit',doraster,issave); %#ok<*ASGLU>
[baseline_fb maxvalue_fb] = poppsth(I,resdir,align2,'hit',doraster,issave); %#ok<*ASGLU>
if issave
    fnm = fullfile(resdir,'stats.mat');
    save(fnm,'baseline_tone','maxvalue_tone','baseline_fb','maxvalue_fb')
end
keyboard

y = log(maxvalue_fb./baseline_fb);  % reliability
x = log(maxvalue_tone./baseline_tone);
[b,bint,r,rint,stats] = regress(y',[ones(length(x),1),x']);
R = sqrt(stats(1));         % correlation coefficient (R-value of the regression)
F = stats(2);           % F-test for H0: all coeff.-s are zero
p = stats(3);           % F-test significance
figure
plot(x,y,'o')
axis square
icp = b(1);   % intercept
gr = b(2);   % gradient
xx = min(x):0.01:max(x);
yy = xx .* gr + icp;
hold on
plot(xx,yy,'Color','k','LineWidth',2)   % overlay regression line
xlabel('Response to tone')
ylabel('Response to water')

% -------------------------------------------------------------------------
function [baseline maxvalue] = poppsth(I,resdir,align,trialtype,doraster,issave)

% Load CellBase if indices to CELLIDLIST are passed
if isnumeric(I)
    loadcb
    I = CELLIDLIST(I);
end

% Time window
wn = [-0.5 0.25];   % in seconds
dt = 0.001;
time = wn(1):dt:wn(2);   % time vector

% Call 'ultimate_psth'
% Preallocate
problem_ids = {};   % cellIDs for which PSTH failed
problem_msg = {};   % corresponding error messages
[baseline maxvalue] = deal([]);  % PSTH stats
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
    switch trialtype
        case 'hit'
            tonelastev = findAlignEvent_posfeedback_gonogo(cellid); %#ok<NASGU>
        case 'fa'
            tonelastev = findAlignEvent_negfeedback_gonogo(cellid); %#ok<NASGU>
    end
    feedbacklastev = []; %#ok<NASGU>
    clastev = eval([align 'lastev']);   % last event is feedback for tone alignment
        
    try
        
        % Calcualte PSTH for Hits
        [psth, spsth, spsth_se, tags, spt, stats] = ultimate_psth(cellid,'trial',cevent,wn,...
            'dt',dt,'display',true,'sigma',0.02,'parts','all',...
            'isadaptive',2,'event_filter','custom','filterinput',cfilter,...
            'first_event',[],'last_event',clastev,...
            'maxtrialno',Inf,'baselinewin',[-0.5 0],'testwin',[0 0.25],...
            'relative_threshold',0.1);
        
        % Concatenate data from different cells
        baseline(k) = stats.baseline;
        maxvalue(k) = stats.maxvalue;
        
        % Save PSTH figure
        H = gcf;
        if issave
            cellidt = regexprep(cellid,'\.','_');
            title([regexprep(cellidt,'_',' ') ' ' align])
            fnm = [resdir cellidt '_PSTH_' align '.fig'];
            saveas(H,fnm)
        end
        close(H)
        
        % Raster plot
        if doraster
            H = figure;
            viewcell2b(cellid,'TriggerName',cevent,'SortEvent',sevent,...
                'eventtype','behav','ShowEvents',{{sevent}},...
                'LastEvents',clastev,'Partitions','#ResponseType','window',wn)
            maximize_figure(H)
            
            % Save raster plot figure
             if issave
                fnm = [resdir cellidt '_raster_' align '.fig'];
                saveas(H,fnm)
                fnm = [resdir cellidt '_raster_' align '.jpg'];
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