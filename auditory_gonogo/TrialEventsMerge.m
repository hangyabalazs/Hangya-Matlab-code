function TE = TrialEventsMerge(datapath,sessionlist)


if nargin < 1,
    %datapath =['/Users/adam/Projects/updating/Rat3'];
    datapath =['/Users/ranades/SoloData/Data/matt/mk05'];
    datapath =['/Users/ranades/Documents/work/Data/behavior/auditory_gonogo/mk05'];
end


list=dir(datapath);
list=setdiff({list.name},{'.';'..';'.DS_Store';'CVS'});

if nargin<2,
    sessionlist = list;
end

% if you pass filenames for select sessions, only those will be used for
% merging.
list = intersect(list,sessionlist);


% list = listfiles(datapath);
TE = [];
% TrialEvents = load(list{1});
TrialEvents = load([datapath filesep list{1}]);

fields =fieldnames(TrialEvents);

ind = 1;

% Index for block number.
bl_ind = 0;
for iL=1:length(list)
    
    if  strcmpi(list{iL}(end-2:end),'mat') && strcmpi(list{iL}(1:2),'TE')
        TrialEvents = load([datapath filesep list{iL}]);
    end
    fields =fieldnames(TrialEvents);
    NumTrials=length(TrialEvents.TrialStart);
    for iF=1:length(fields)
        if length(TrialEvents.(fields{iF})) == NumTrials
            TE.(fields{iF})(ind:ind+NumTrials-1) = TrialEvents.(fields{iF})(1:NumTrials);
        end
        TE.SessionID(ind:ind+NumTrials-1) = iL;
    end %iF
    TE.BlockNum(ind:ind+NumTrials -1) = TE.BlockNum(ind:ind+NumTrials -1)+bl_ind;
    ind = ind+NumTrials;
    bl_ind = max(TE.BlockNum);
end %iL