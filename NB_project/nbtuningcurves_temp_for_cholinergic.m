function nbtuningcurves(I,issave)
%NBTUNINGCURVES   PSTHs for different tone intensitied.
%   NBTUNINGCURVES(I,ISSAVE) calculates non-adaptive PSTHs for a set of
%   cells aligned to 'Go' and 'No-go' response onset ('Hit' and 'False
%   Alarm'). The trials are partitioned based on stimulus intensity. The
%   PSTHs are plotted and saves.
%   Input parameters: 
%       I - index set to CELLIDLIST (see CellBase documentation); if empty
%           or not specified, a predefined set of basal forebrain cells is
%           used.
%       ISSAVE - controls saving
%
%   See also ULTIMATE_PSTH and NBRESPONSESORTER.

%   Edit log: BH 9/5/12
 
% Pass the control to the user in case of error
dbstop if error
choosecb('NB')

% Input argument check
error(nargchk(0,2,nargin))
if nargin < 2
    issave = true;
end
if nargin < 1
    I = [];
end
 
% Directories
global DATAPATH
fs = filesep;
resdir = [DATAPATH 'NB' fs 'tuningcurves3_PV' fs];
 
% List of cellIDs
if isempty(I)
    selstr = ['"ChAT+"==1&"ID_PC">20&"Lr_PC"<0.15&' ...
        'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})'];
    ChAT = selectcell(selstr);
    I = cellfun(@findcellpos,ChAT);
    
    pChAT = {'n029_120210a_3.3' 'n029_120214b_2.7' 'n029_120215a_2.2' ...
        'n029_120215a_3.4' 'n029_120221b_6.1' 'n029_120222b_4.1'};
    I = cellfun(@findcellpos,pChAT);
    
    selstr = ['"PV+"==1&"ID_PC">20&"Lr_PC"<0.15&' ...
        'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})'];
    PV = selectcell(selstr);
    I = cellfun(@findcellpos,PV);
end
I = I(:)';   % convert to row vector

% PSTH
poppsth(I,resdir,issave); %#ok<*ASGLU>

% -------------------------------------------------------------------------
function poppsth(I,resdir,issave)

% Load CellBase
loadcb

% Time window
wn = [-0.6 0.6];   % in seconds
dt = 0.001;
time = wn(1):dt:wn(2);   % time vector

% Call 'ultimate_psth'
% Preallocate
problem_ids = {};
problem_msg = {};
for k = I
    cellid = CELLIDLIST{k}; %#ok<USENS>
    disp(cellid)
    try
        
        % Calcualte PSTH for Hits
        [psth, spsth, spsth_se, tags] = ultimate_psth(cellid,'trial','LeftWaterValveOn',wn,...
            'dt',dt,'display',true,'sigma',0.005,'parts','#StimulusDuration',...
            'isadaptive',0,'event_filter','custom','filterinput','Hit==1',...
            'maxtrialno',Inf,'baselinewin',[-0.5 0],'testwin',[0 0.5],...
            'relative_threshold',0.1);
        
        % Plot
        figure
        hold on
        NumPsth = size(spsth,1);
        clr = summer(NumPsth);
        for pk = 1:NumPsth
%             plot(spsth(pk,:),'Color',clr(pk,:),'LineWidth',2)
            errorshade(time,spsth(pk,:),spsth_se(pk,:),'LineWidth',2,'LineColor',clr(pk,:),...
                'ShadeColor',clr(pk,:))
        end
        legend(tags)
        xlim([-0.2 0.2])
        
        % Save figure
        if issave
            H = gcf;
            cellidt = regexprep(cellid,'\.','_');
            title([regexprep(cellidt,'_',' ') ' Hit'])
            fnm = [resdir cellidt '_PSTH_Hit.fig'];
            saveas(H,fnm)
        end
        close(H)
        
        % Calcualte PSTH for FAs
        [psth, spsth, ~, tags] = ultimate_psth(cellid,'trial','LeftPortIn',wn,...
            'dt',0.001,'display',true,'sigma',0.005,'parts','#StimulusDuration',...
            'isadaptive',0,'event_filter','custom','filterinput','FalseAlarm==1',...
            'maxtrialno',Inf,'baselinewin',[-0.5 0],'testwin',[0 0.5],...
            'relative_threshold',0.1);
        
        % Plot
        figure
        hold on
        NumPsth = size(spsth,1);
        clr = summer(NumPsth);
        for pk = 1:NumPsth
%             plot(spsth(pk,:),'Color',clr(pk,:),'LineWidth',2)
            errorshade(time,spsth(pk,:),spsth_se(pk,:),'LineWidth',2,'LineColor',clr(pk,:),...
                'ShadeColor',clr(pk,:))
        end
        legend(tags)
        xlim([-0.2 0.2])
        
        % Save figure
        if issave
            H = gcf;
            cellidt = regexprep(cellid,'\.','_');
            title([regexprep(cellidt,'_',' ') ' False Alarm'])
            fnm = [resdir cellidt '_PSTH_FA.fig'];
            saveas(H,fnm)
        end
        close(H)
        
    catch ME
        
        % Error handling
        problem_ids{end+1} = cellid;  %#ok<AGROW> % collect cellIDs resulting in error
        problem_msg{end+1} = ME.message;  %#ok<AGROW> % store corresponding error messages
        disp(['Something went wrong for cell ' cellid '.'])
        disp(ME.message)
    end
end