function bimodalITI7
%BIMODALITI7   Reaction time analysis.
%   BIMODALITI7 analyzes reaction (RT) time in the function of the
%   foreperiod. All sessions with bimodal foreperiod distributions are
%   included; trials are pooled from tese sessions. BIMODALITI7 plots
%   median RT against foreperiod length, separated according to stimulus
%   intensity and smoothed with a moving average.
%
%   See also BIMODALITI6.

% Specify folder
basedir = getpref('cellbase','datapath');
mousename = 'n023';

% Sessions
% Animal 'n023'
SESSIONS.n023.bimodal0 = {'111220a' '111221a' '111222a' '111223a' '111227a' '111228a' ...
    };    % bimodal, different mixing proportions (narrow peaks)
SESSIONS.n023.bimodal = {'111229a' '111230a' '111231a' '120101a' '120102a' '120103a' ...
    '120104a' '120105a' '120106a' '120110a' '120111a' '120112a' '120113a' ...
    '120114a' '120115a' '120116a'};    % bimodal
SESSIONS.n023.unimodal0 = {'120117a' '120118a' '120119a' '120120a' '120121a' '120122a'};   % unimodal
SESSIONS.n023.unimodal = {'120123a' '120124a' '120125a' '120126a' '120127a' '120128a'...
    '120131a' '120201a' '120202a' '120203a'};   % unimodal with uniform
SESSIONS.n023.exponential = {'120207a' '120208a' '120209a' '120210a' '120211a' ...
    '120213a' '120214a' '120215a' '120216a' '120217a'};     % exponential
SESSIONS.n023.unimodalnew = {'120220a' '120221a'};       % unimodal again

% Animal 'n026'
SESSIONS.n026.exponential = {'120207a' '120208a' '120209a' '120210a' '120211b' ...
    '120213a' '120214a' '120215a' '120216a'};     % exponential
SESSIONS.n026.bimodal = {'120218a' '120220a' '120221a' '120222a' '120223a'};  % bimodal (not performing after Cosyne, not included)
        
% Animal 'n027'
SESSIONS.n027.exponential = {'120208a' '120209a' '120210a' '120211a'};     % exponential
SESSIONS.n027.bimodal = {'120213b' '120214a' '120215a' '120215c' '120216a' '120217a' ...
    '120218a' '120220a' '120221a' '120222a' '120223a' ...
    '120301a' '120302a' '120303a' '120305a' '120306a'...
    '120307a' '120308a' '120309a' '120310a' '120311a'...
    '120313a' '120314a' '120315a'};     % bimodal (120215a and c: not doing)
        
% Animal 'n028'
SESSIONS.n028.exponential = {'120207a' '120208a' '120210a' '120211a' ...
    '120213a' '120214a' '120215a' '120216a'};     % exponential
SESSIONS.n028.bimodal = {'120217a' '120218a' '120220a' '120221a'...
    '120222a' '120222c' '120223a' ...
    '120301a' '120302a' '120303a' '120305a' '120306a' '120307a'};     % bimodal
        
% Animal 'n029'
SESSIONS.n029.bimodal = {'120204a' '120204b' '120206a' '120207a' '120207b' ...
    '120208a' '120209a' '120210a' '120211a' '120213a'...
    '120213b' '120214a' '120214b' '120215a' '120216b' '120217a' ...
    '120218a' '120220a' '120220b' '120221a' '120221b' '120222a' '120222b' ...
    '120301a' '120302a' '120303a' '120305a'...
    '120306a' '120307a' '120308a' '120309a' '120310a' '120311a'...
    '120313a' '120314a' '120315a'};     % bimodal (remark: 120209a session very long and good RT data!!!)

% Animal 'n037'
SESSIONS.n037.bimodal0 = {'120909a' '120910a' '120911a' '120911b' ...  % before STIM session, not too good
    '120913a' '120915a' '120916a' '120917a' '120918a' '120919a'};  % bimodal with narrow peaks
SESSIONS.n037.bimodal = {'120920b' '120921a' '120922a' '120923a' '120923d' ...
    '120924a' '120925a' '120926a' '120928a' '121001a' '121002a' '121003a' ...
    '121004a' '121005b' '121006a' '121008a' '121009a' '121010a' '121011a' ...
    '121012a'};     % bimodal

% Animal 'n038'
SESSIONS.n038.bimodal0 = {'120909a' '120910a' '120911a' '120912a' ...
    '120913b' '120914a' '120917a' '120918a' '120919a'};  % bimodal with narrow peaks
SESSIONS.n038.bimodal = {'120920a' '120921a' '120922a' ...
    '120924a' '120925a' '120926a' '120928a' '121001a' '121002b' '121003a' ...
    '121004a' '121005a' '121006a' '121008a' '121009a' '121010a' '121011a' ...
    '121012a'};     % bimodal

% Animal 'n039'
SESSIONS.n039.bimodal0 = {'120909a' '120910a' '120911a' '120912a' ...
    '120913a' '120914a' '120917b' '120918a' '120919a' '120919b' '120920a'};  % bimodal with narrow peaks
SESSIONS.n039.bimodal = {'120921b' ...
    '120924a' '120925a' '120926a' '120928a' '120929a' '121001a' ...
    '121004a'};     % bimodal

% Animal 'n040'
SESSIONS.n040.bimodal0 = {'121010a' '121011a' '121012a' '121015a' ...
    '121016a' '121019a' '121021a' '121022a'};  % bimodal with narrow peaks
SESSIONS.n040.bimodal = {'121023a' '121024a' '121025a' '121026a' '121026c' '121027a' ...
    '121108a' '121110a' '121111a' '121112a' '121113a' '121114a' '121115a' ...
    '121119a' '121120a'};     % bimodal

% Animal 'n043'
SESSIONS.n043.bimodal0 = {'121125a' '121126a' '121127a' ...
    '121204a' '121205a' '121206a' '121207a' '121210a' '121211a' ...
    '121213a' '121214a'};  % bimodal with narrow peaks
SESSIONS.n043.bimodal = {'121217a' '121218a' '121219a' '121220a' '121221a'};     % bimodal

% Animal 'n045'
SESSIONS.n045.bimodal0 = {'121210a' '121211a' '121212a' '121213a' '121214a' '121217a'};  % bimodal with narrow peaks
SESSIONS.n045.bimodal = {'121217b' '121218a' '121219a' '121220a' '121221a' ...
    '121229a' '121231a' '130102a' '130104a' '130107a' '130109a'};     % bimodal

% Animal 'n046'
SESSIONS.n046.bimodal0 = {'121207a' '121210a' '121211a' '121212a' '121213a' '121214a' ...
    '121217a' '121218a'};  % bimodal with narrow peaks
SESSIONS.n046.bimodal = {'121219a' '121220a' '121221a' ...
    '121229a' '121229b' '121229c' '121229d' '121229g' '121230a' ...
    '121231a' '130101a' '130102a' '130102b' '130102d' '130102e' ...
    '130103a' '130104a' '130107a' '130108a' '130108d' '130108e' '130108f' ...
    '130109a' '130110a' '130111a'};     % bimodal

% Subject 'mf'
SESSIONS.mf.bimodal = {'120806a' '120806b' '120806c' '120807a' ...
    '120807b' '120807c' '120808a' '120808b' '120808c' '120809a'};
SESSIONS.mf.bimodal2 = {'120814a' '120823a' '120827a'};

% Subject 'ak'
SESSIONS.ak.bimodal = {'120815a' '120815b' '120815c' '120816a' '120816b' ...
    '120816c' '120817a' '120821a' '120821b' '120821c' '120822a' '120822b'};

% Subject 'ld'
SESSIONS.ld.bimodal = {'121115b' '121126a' '121210a' '121218a' '121228a' ...
    '130114a' '130203a' '130211a'};    %#ok<STRNU>

eval(['sessions = SESSIONS.' mousename '.exponential;']);
numSessions = length(sessions); %#ok<USENS>

% Preallocate
meanrt = cell(numSessions,1);
meanhit = cell(numSessions,1);
meanmiss = cell(numSessions,1);
medianrt = cell(numSessions,1);
medianhit = cell(numSessions,1);
medianmiss = cell(numSessions,1);
nrt = cell(numSessions,1);
nhit = cell(numSessions,1);
nhits = cell(numSessions,1);
nmiss = cell(numSessions,1);
nmisss = cell(numSessions,1);
allrt = cell(numSessions,1);
alliti = cell(numSessions,1);
allhit = cell(numSessions,1);
allmiss = cell(numSessions,1);
allsi = cell(numSessions,1);

allrt2 = [];
alliti2 = [];
allhit2 = [];
allmiss2 = [];
allsi2 = [];

% Binning
edges = 0.1:0.25:3;   % mouse bimodal
% edges = 0.1:0.5:7.5;  % for human bimodal, longer times
cnts = (edges(1:end-1) + edges(2:end)) / 2;
bno = length(cnts);

% Get session data
ints = [20 30 40 50];
if isequal(whichcb,'human_susattn')
    ints = [-3.2189 -2.6912 -2.1637 -1.6368];   % intensities for human subject 'mf'
end
if isequal(mousename,'ld')
    ints = [-3.4762 -3.2130 -2.9499 -2.6867];   % intensities for human subject 'ld'
end
lenints = length(ints);
clr = summer(lenints);
H1 = figure;
for SoundInt = 1:lenints
    for iS = 1:numSessions
        session = sessions{iS};
        datapath = fullfile(basedir,mousename,session);
        list = dir(datapath);
        list = setdiff({list.name},{'.';'..';'.DS_Store'});
        stc = strfind(list,'TE.mat');
        cf = ~cellfun(@isempty,stc);
        list = list(cf);
        
        % Load behavioral data
        TE = load(fullfile(datapath,list{1}));
        
        % Maintain compatibility between human and mouse version of the task
        if isequal(whichcb,'human_susattn')
            TE.StimulusDuration = log(TE.SoundIntensity);
            TE.StimulusDuration = ...
                round(TE.StimulusDuration*10000) / 10000;   % limit precision of tone intensities to 4 digits
        end
    
    % Find unique sound intensities
        SI = unique(TE.StimulusDuration);
        SIn = SI(~isnan(SI));
        if ~isequal(SIn,ints)
            warning('bimodalITI:differentIntensities','Intensities used differ from default.')
        end
        
        % RT vs ITI
        TE2 = filterTE(TE,TE.StimulusDuration==ints(SoundInt));
        [meanrt{iS} medianrt{iS} nrt{iS} allrt{iS} alliti{iS}] = condmean([TE2.ITIDistribution' TE2.GoRT'],edges);
        [meanhit{iS} medianhit{iS} nhit{iS} allhit{iS} alliti{iS}] = condmean([TE2.ITIDistribution' TE2.Hit'],edges);
        nhits{iS} = cellfun(@(g)sum(nan2zero(g)),allhit{iS});
        [meanmiss{iS} medianmiss{iS} nmiss{iS} allmiss{iS} alliti{iS}] = condmean([TE2.ITIDistribution' TE2.Miss'],edges);
        nmisss{iS} = cellfun(@(g)sum(nan2zero(g)),allmiss{iS});
        [~, ~, ~, allsi{iS}] = condmean([TE2.ITIDistribution' TE2.StimulusDuration'],edges);
        allrt2 = [allrt2 TE2.GoRT];
        alliti2 = [alliti2 TE2.ITIDistribution];
        allhit2 = [allhit2 TE2.Hit];
        allmiss2 = [allmiss2 TE2.Miss];
        allsi2 = [allsi2 TE2.StimulusDuration];
    end
    
    % Pool all data bin-wise and calculate median
    mnsall = cell(1,bno);
    mnsmedian = nan(1,bno);
    mnsmedian_se = nan(1,bno);
    for t = 1:bno
        for k = 1:numSessions
            mnsall{t} = [mnsall{t}; allrt{k}{t}];
        end
        mnsmedian(t) = nanmedian(mnsall{t});
        if ~isempty(mnsall{t}) && ~all(isnan(mnsall{t}))
            mnsmedian_se(t) = se_of_median(mnsall{t},1000);
        end
    end
    
    % Plot smoothed version
    smkernel = 3;
%     smkernel = 1;
    smnsmedian = smooth(mnsmedian,'linear',smkernel);
%     smnsmedian_se = smooth(mnsmedian_se,'linear',smkernel);
    figure(H1);
    plot(cnts,smnsmedian,'Color',clr(SoundInt,:),'LineWidth',2)
    hold on
%     errorshade(cnts,smnsmedian,smnsmedian_se,...
%         'LineColor',clr(SoundInt,:),'ShadeColor',clr(SoundInt,:))
    [nm xout] = hist(TE.ITIDistribution,50);
    stairs(xout,nm/900+0.18,'Color','r')
    
%     keyboard
end

% -------------------------------------------------------------------------
function TE = filterTE(TE,valid_trials)

fnm = fieldnames(TE);
for k = 1:length(fnm)
    TE.(fnm{k}) = TE.(fnm{k})(valid_trials);
end