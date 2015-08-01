function training_vs_hitresponse
%TRAINING_VS_HITRESPONSE   Correlation between training history and reward-related firing rate increase.
%   TRAINING_VS_HITRESPONSE plots relative firing rate change, i.e. the 
%   logarithm of maximal firing rate divided by baseline firing rate as a 
%   function of training parameters for all cholinergic and putative
%   cholinergic neurons (see NBCLUSTERING). Correlation coefficient and
%   significance of regression are calculated.
%
%   See also DEPTH_VS_HITRESPONSE, NBRESPONSESORTER and NBCLUSTERING.

% Cholinergic neurons
ChAT = selectcell(['"ChAT+"==1&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("session_type",{''gonogo'',''feedbackdelay''})&' ...
    'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})']);  % identified
ChAT = setdiff(ChAT,'n072_141222a_4.3');   % hit and FA resonse cannot be compared because of strong drift
ChAT = [ChAT 'n045_121217x_4.6'];   % clustered based on light-evoked spikes
pChAT = selectcell(['"pChAT+"==1&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("session_type",{''gonogo'',''feedbackdelay''})&' ...
    'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})']);  % putative
% allChAT = [ChAT pChAT];   % identified and putative
allChAT = ChAT;
NumChAT = length(allChAT);   % number of cholinergic cells

% Color-marker code
mouse = getvalue('RatID',allChAT);  % animal ID
nms = length(unique(mouse));   % number of mice (n=7)
mouse2 = nan(size(mouse));
for k = 1:nms
    mouse2(mouse==min(mouse)) = k;
    mouse(mouse==min(mouse)) = Inf;
end
clr = colormap(jet(nms));   % unique color for each mouse
mrk = 'Sp^hdo<>VX*+';   % unique marker for each mouse
colorcode = false;   % control color code
markercode = true;   % control marker code

% Training history
[hitnum fanum trialnum sessionnum daynum] = deal(nan(NumChAT,1));
for iC = 1:NumChAT
    training_history = trhist(allChAT{iC});
    nums = [training_history.Hits];
    hitnum(iC) = sum(nums);   % number of hits in prev. sessions
    nums = [training_history.FAs];
    fanum(iC) = sum(nums);   % number of false alarms in prev. sessions
    nums = [training_history.Trials];
    trialnum(iC) = sum(nums);   % number of trials in prev. sessions
    sessionnum(iC) = length(training_history);   % number of prev. sessions
    daynum(iC) = training_history(end).Days;   % number of days since first session
end

% Reward related firing rate change
Hitstat = getvalue('Hit_psth_stats',allChAT);   % hit trial PSTH statistics
HitPSTH = getvalue('Hit_psth',allChAT);   % hit trial PSTHs
[Hiteff Hiteff2 baselineH] = deal(nan(1,NumChAT));
FAstat = getvalue('FA_psth_stats',allChAT);   % hit trial PSTH statistics
FAPSTH = getvalue('FA_psth',allChAT);   % hit trial PSTHs
[FAeff FAeff2 baselineF] = deal(nan(1,NumChAT));
RelEff = nan(1,NumChAT);

xvar = daynum;
figure
hold on
for k = 1:NumChAT
%     disp(allChAT{k})
%     pause(2)
    Hiteff(k) = Hitstat{k}.maxvalue;   % maximal FR
    baselineH(k) = Hitstat{k}.baseline;  % baseline FR
    Hiteff2(k) = log(Hiteff(k)/baselineH(k));  % relative FR change on log scale
    FAeff(k) = FAstat{k}.maxvalue;   % maximal FR
    baselineF(k) = FAstat{k}.baseline;  % baseline FR
    FAeff2(k) = log(FAeff(k)/baselineF(k));  % relative FR change on log scale
    if colorcode
        cclr = clr(mouse2(k),:);   % unique color for each mouse
    else
        lc = length(ChAT);   % number of identified neurons
        cclr = [0 0.8*(k<=lc) 0.8*(k>lc)];   % color according to identified/putative
    end
    if markercode
        cmrk = mrk(mouse2(k));   % unique marker for each mouse
    else
        cmrk = 'o';
    end
    RelEff(k) = log(max([Hiteff(k)-baselineH(k) 0.5]))-log(FAeff(k)-baselineF(k));
    plot(xvar(k),RelEff(k),cmrk,'MarkerSize',12,'MarkerFaceColor',cclr,'MarkerEdgeColor','k')  % scatter plot
end

% Regression
y = RelEff;
x = xvar';
[b,bint,r,rint,stats] = regress(y',[ones(length(x),1),x']);
R = sqrt(stats(1));         % correlation coefficient (R-value of the regression)
F = stats(2);           % F-test for H0: all coeff.-s are zero
p = stats(3);           % F-test significance
xlabel('Depth (um)')
setmyplot_balazs
axis square
icp = b(1);   % intercept
gr = b(2);   % gradient
xx = min(x):0.01:max(x);
yy = xx .* gr + icp;
hold on
plot(xx,yy,'Color','k','LineWidth',2)   % overlay regression line
text('Units','normalized','Position',[0.7 0.7],...
    'String',{['p = ' num2str(p)] ['R = ' num2str(R)]})

keyboard

% -------------------------------------------------------------------------
function training_history = trhist(cellid)

% Sessions
[mouse session] = cellid2tags(cellid);   % animal ID
mousedir = fullfile(getpref('cellbase','datapath'),mouse);   % corresponding directory in CellBase
dmd = dir(mousedir);
dmd = dmd(3:end);
sessions = {dmd.name};   % session names

% Session info
NumSessions = length(sessions);   % number of sessions
next = 0;
training_history = struct('Trials',{},'Hits',{},'FAs',{},'SessionID',{},'Days',{});
for iS = 1:NumSessions
    if isequal(sessions{iS},session)
        break   % we only need the sessions BEFORE the cell was recorded
    end
    sessionpath = fullfile(mousedir,sessions{iS});   % session pathname
    behav_file = [sessionpath filesep 'TE.mat'];   % trial events file name
    if exist(behav_file,'file')
        TE = load(behav_file);   % trial events structure
        if isfield(TE,'FalseAlarm') && nansum(TE.FalseAlarm) > 0 && ...  % session with full behav. task
                ismember(sessionpath(end),...
                {'a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p'})   % other letters may correspond to re-clustered sessions
            next = next + 1;
            if next == 1
                first_session = sessions{iS};   % first session with full task
            end
            training_history(next).Trials = length(TE.Hit);   % number of trials
            training_history(next).Hits = nansum(TE.Hit);   % number of hits
            training_history(next).FAs = nansum(TE.FalseAlarm);   % number of false alarms
            training_history(next).SessionID = sessions{iS};   % session ID
            
            session_date = datenum(str2double(['20' session(1:2)]),str2double(session(3:4)),str2double(session(5:6)));  % date of the recording session
            first_session_date = datenum(str2double(['20' first_session(1:2)]),str2double(first_session(3:4)),str2double(first_session(5:6)));   % date of the first session with full task
            training_history(next).Days = session_date - first_session_date;  % days since first session
        end
    end
end