%% first attempt

% Load
global DATADIR
global DATAPATH
load([DATADIR 'SOM_Sachin\PSTH\tagged_list_agreed']);
loadcb
load([DATAPATH 'SOM\ResponseCCG3\CCG_STRENGTH.mat'])

% Subsets
[is som_inx inx2] = intersect(CELLIDLIST,som_beh);
[is pv_inx inx2] = intersect(CELLIDLIST,pv_beh);
nt_inx = setdiff((1:length(CELLIDLIST)),[som_inx pv_inx]);

% Plot summary
UaStrength = squeeze(CcgStrength(som_inx,:,4));
mx = max(UaStrength,[],2);
for k = 1:length(som_inx)
    if ~isequal(sum(UaStrength(k,:)),0)
        UaStrength(k,:) = UaStrength(k,:)==mx(k);
    end
end
figure;
plot(sum(UaStrength)/sum(UaStrength(:)),'Color','b')

UaStrength = squeeze(CcgStrength(pv_inx,:,4));
mx = max(UaStrength,[],2);
for k = 1:length(pv_inx)
    if ~isequal(sum(UaStrength(k,:)),0)
        UaStrength(k,:) = UaStrength(k,:)==mx(k);
    end
end
hold on
plot(sum(UaStrength)/sum(UaStrength(:)),'Color','r')

UaStrength = squeeze(CcgStrength(nt_inx,:,4));
mx = max(UaStrength,[],2);
for k = length(nt_inx)
    if ~isequal(sum(UaStrength(k,:)),0)
        UaStrength(k,:) = UaStrength(k,:)==mx(k);
    end
end
plot(sum(UaStrength)/sum(UaStrength(:)),'Color',[0.7 0.7 0.7])

%% second attempt

% Load
global DATADIR
global DATAPATH
load([DATADIR 'SOM_Sachin\PSTH\tagged_list_agreed']);
loadcb
load([DATAPATH 'SOM\ResponseCCG3\CCG_STRENGTH.mat'])

% Subsets
[is som_inx inx2] = intersect(CELLIDLIST,som_beh);
[is pv_inx inx2] = intersect(CELLIDLIST,pv_beh);
nt_inx = setdiff((1:length(CELLIDLIST)),[som_inx pv_inx]);

% Plot summary
UaS = squeeze(CcgStrength(som_inx,:,2));
LaS = squeeze(CcgStrength(som_inx,:,3));
UaStrength = squeeze(CcgStrength(som_inx,:,4));
for k = 1:length(som_inx)
    tuas = UaStrength(k,:);
    tu = UaS(k,:);
    tl = LaS(k,:);
    tinx = [];
    if ~isequal(sum(tuas),0)
        tinx = find(tu==0&tl==0);
        tuas(tinx) = zeros(1,length(tinx));
        UaStrength(k,:) = tuas;
    end
end
figure;
plot(sum(UaStrength)/sum(UaStrength(:)),'Color','b')

UaStrength = squeeze(CcgStrength(pv_inx,:,4));
mx = max(UaStrength,[],2);
for k = 1:length(pv_inx)
    if ~isequal(sum(UaStrength(k,:)),0)
        UaStrength(k,:) = UaStrength(k,:)==mx(k);
    end
end
hold on
plot(sum(UaStrength)/sum(UaStrength(:)),'Color','r')

UaStrength = squeeze(CcgStrength(nt_inx,:,4));
mx = max(UaStrength,[],2);
for k = length(nt_inx)
    if ~isequal(sum(UaStrength(k,:)),0)
        UaStrength(k,:) = UaStrength(k,:)==mx(k);
    end
end
plot(nansum(UaStrength)/nansum(UaStrength(:)),'Color',[0.7 0.7 0.7])


%% third attempt - max

% Load
global DATADIR
global DATAPATH
load([DATADIR 'SOM_Sachin\PSTH\tagged_list_agreed']);
loadcb
load([DATAPATH 'SOM\ResponseCCG3\CCG_STRENGTH.mat'])

% Subsets
[is som_inx inx2] = intersect(CELLIDLIST,som_beh);
[is pv_inx inx2] = intersect(CELLIDLIST,pv_beh);
nt_inx = setdiff((1:length(CELLIDLIST)),[som_inx pv_inx]);

% Plot summary
UaStrength = squeeze(CcgStrength(som_inx,:,4));
mx = max(UaStrength,[],2);
for k = 1:length(som_inx)
    UaStrength(k,:) = UaStrength(k,:)==mx(k);
end
figure;
plot(sum(UaStrength)/sum(UaStrength(:)),'Color','b')

UaStrength = squeeze(CcgStrength(pv_inx,:,4));
mx = max(UaStrength,[],2);
for k = 1:length(pv_inx)
    UaStrength(k,:) = UaStrength(k,:)==mx(k);
end
hold on
plot(sum(UaStrength)/sum(UaStrength(:)),'Color','r')

UaStrength = squeeze(CcgStrength(nt_inx,:,4));
mx = max(UaStrength,[],2);
for k = 1:length(nt_inx)
    UaStrength(k,:) = UaStrength(k,:)==mx(k);
end
plot(nansum(UaStrength)/nansum(UaStrength(:)),'Color',[0.7 0.7 0.7])

%% 4th attempt - rank

% Load
global DATADIR
global DATAPATH
load([DATADIR 'SOM_Sachin\PSTH\tagged_list_agreed']);
loadcb
load([DATAPATH 'SOM\ResponseCCG3\CCG_STRENGTH.mat'])

% Subsets
[is som_inx inx2] = intersect(CELLIDLIST,som_beh);
[is pv_inx inx2] = intersect(CELLIDLIST,pv_beh);
nt_inx = setdiff((1:length(CELLIDLIST)),[som_inx pv_inx]);

% Plot summary
UaStrength = squeeze(CcgStrength(som_inx,:,4));
for k = 1:length(som_inx)
    UaStrength(k,:) = tiedrank(UaStrength(k,:))-1;
end
figure;
plot(sum(UaStrength)/sum(UaStrength(:)),'Color','b')

UaStrength = squeeze(CcgStrength(pv_inx,:,4));
for k = 1:length(pv_inx)
    UaStrength(k,:) = tiedrank(UaStrength(k,:))-1;
end
hold on
plot(sum(UaStrength)/sum(UaStrength(:)),'Color','r')

UaStrength = squeeze(CcgStrength(nt_inx,:,4));
for k = 1:length(nt_inx)
    UaStrength(k,:) = tiedrank(UaStrength(k,:))-1;
end
plot(nansum(UaStrength)/nansum(UaStrength(:)),'Color',[0.7 0.7 0.7])


%% 5th attempt - mean SD

% Load
global DATADIR
global DATAPATH
load([DATADIR 'SOM_Sachin\PSTH\tagged_list_agreed']);
loadcb
load([DATAPATH 'SOM\ResponseCCG3\CCG_STRENGTH.mat'])

% Subsets
[is som_inx inx2] = intersect(CELLIDLIST,som_beh);
[is pv_inx inx2] = intersect(CELLIDLIST,pv_beh);
nt_inx = setdiff((1:length(CELLIDLIST)),[som_inx pv_inx]);

% Plot summary
UaStrength = squeeze(CcgStrength(som_inx,:,4));
figure;
plot(mean(UaStrength),'Color','b')

UaStrength = squeeze(CcgStrength(pv_inx,:,4));
hold on
plot(mean(UaStrength),'Color','r')

UaStrength = squeeze(CcgStrength(nt_inx,:,4));
plot(nanmean(UaStrength),'Color',[0.7 0.7 0.7])

%% 6th attempt - median SD

% Load
global DATADIR
global DATAPATH
load([DATADIR 'SOM_Sachin\PSTH\tagged_list_agreed']);
loadcb
load([DATAPATH 'SOM\ResponseCCG3\CCG_STRENGTH.mat'])

% Subsets
[is som_inx inx2] = intersect(CELLIDLIST,som_beh);
[is pv_inx inx2] = intersect(CELLIDLIST,pv_beh);
nt_inx = setdiff((1:length(CELLIDLIST)),[som_inx pv_inx]);

% Plot summary
UaStrength = squeeze(CcgStrength(som_inx,:,4));
figure;
plot(median(UaStrength),'Color','b')

UaStrength = squeeze(CcgStrength(pv_inx,:,4));
hold on
plot(median(UaStrength),'Color','r')

UaStrength = squeeze(CcgStrength(nt_inx,:,4));
plot(nanmedian(UaStrength),'Color',[0.7 0.7 0.7])

%% 7th attempt - mean SD, normalized

% Load
global DATADIR
global DATAPATH
load([DATADIR 'SOM_Sachin\PSTH\tagged_list_agreed']);
loadcb
load([DATAPATH 'SOM\ResponseCCG3\CCG_STRENGTH.mat'])

% Subsets
[is som_inx inx2] = intersect(CELLIDLIST,som_beh);
[is pv_inx inx2] = intersect(CELLIDLIST,pv_beh);
nt_inx = setdiff((1:length(CELLIDLIST)),[som_inx pv_inx]);

% Plot summary
UaStrength = squeeze(CcgStrength(som_inx,:,4));
SumCcg = SumCcgs(som_inx);
for k = 1:length(som_inx)
    UaStrength(k,:) = (UaStrength(k,:)) / SumCcg(k);
end
figure;
plot(mean(UaStrength),'Color','b')

UaStrength = squeeze(CcgStrength(pv_inx,:,4));
SumCcg = SumCcgs(pv_inx);
for k = 1:length(pv_inx)
    UaStrength(k,:) = (UaStrength(k,:)) / SumCcg(k);
end
hold on
plot(mean(UaStrength),'Color','r')

UaStrength = squeeze(CcgStrength(nt_inx,:,4));
SumCcg = SumCcgs(nt_inx);
for k = 1:length(nt_inx)
    UaStrength(k,:) = (UaStrength(k,:)) / SumCcg(k);
end
plot(nanmean(UaStrength),'Color',[0.7 0.7 0.7])

%% 8th attempt - mean SD, normalized

% Load
global DATADIR
global DATAPATH
load([DATADIR 'SOM_Sachin\PSTH\tagged_list_agreed']);
loadcb
load([DATAPATH 'SOM\ResponseCCG3\CCG_STRENGTH.mat'])

% Subsets
[is som_inx inx2] = intersect(CELLIDLIST,som_beh);
[is pv_inx inx2] = intersect(CELLIDLIST,pv_beh);
nt_inx = setdiff((1:length(CELLIDLIST)),[som_inx pv_inx]);

% Plot summary
UaStrength = squeeze(CcgStrength(som_inx,:,4));
for k = 1:length(som_inx)
    UaStrength(k,:) = (UaStrength(k,:)) / sum(UaStrength(k,:));
end
figure;
plot(mean(UaStrength),'Color','b')

UaStrength = squeeze(CcgStrength(pv_inx,:,4));
SumCcg = SumCcgs(pv_inx);
for k = 1:length(pv_inx)
    UaStrength(k,:) = (UaStrength(k,:)) / sum(UaStrength(k,:));
end
hold on
plot(mean(UaStrength),'Color','r')

UaStrength = squeeze(CcgStrength(nt_inx,:,4));
SumCcg = SumCcgs(nt_inx);
for k = 1:length(nt_inx)
    UaStrength(k,:) = (UaStrength(k,:)) / sum(UaStrength(k,:));
end
plot(nanmean(UaStrength),'Color',[0.7 0.7 0.7])

%% 9th attempt - max for sign.

% Load
global DATADIR
global DATAPATH
load([DATADIR 'SOM_Sachin\PSTH\tagged_list_agreed']);
loadcb
load([DATAPATH 'SOM\ResponseCCG3\CCG_STRENGTH.mat'])

% Subsets
[is som_inx inx2] = intersect(CELLIDLIST,som_beh);
[is pv_inx inx2] = intersect(CELLIDLIST,pv_beh);
nt_inx = setdiff((1:length(CELLIDLIST)),[som_inx pv_inx]);

% Plot summary
UaS = squeeze(CcgStrength(:,:,2));
LaS = squeeze(CcgStrength(:,:,3));
UaStrength = squeeze(CcgStrength(:,:,4));
for k = 1:size(CcgStrength,1)
    tuas = UaStrength(k,:);
    tu = UaS(k,:);
    tl = LaS(k,:);
    tinx = tu<0.0001&tl<0.0001;
    tuas(tinx) = NaN;
    ic1 = tuas==nanmax(tuas);
    ic2 = ~isnan(tuas) & ~ic1;
    tuas(ic1) = 1;
    tuas(ic2) = 0;
    UaStrength(k,:) = tuas;
end

UaStrength1 = UaStrength(som_inx,:);
figure;
plot(nansum(UaStrength1)/nansum(UaStrength1(:)),'Color','b')

UaStrength1 = UaStrength(pv_inx,:);
hold on
plot(nansum(UaStrength1)/nansum(UaStrength1(:)),'Color','r')

UaStrength1 = UaStrength(nt_inx,:);
plot(nansum(UaStrength1)/nansum(UaStrength1(:)),'Color',[0.7 0.7 0.7])

%% 10th attempt - mean norm. SD for sign.

% Load
global DATADIR
global DATAPATH
load([DATADIR 'SOM_Sachin\PSTH\tagged_list_agreed']);
loadcb
load([DATAPATH 'SOM\ResponseCCG3\CCG_STRENGTH.mat'])

% Subsets
[is som_inx inx2] = intersect(CELLIDLIST,som_beh);
[is pv_inx inx2] = intersect(CELLIDLIST,pv_beh);
nt_inx = setdiff((1:length(CELLIDLIST)),[som_inx pv_inx]);

% Plot summary
UaS = squeeze(CcgStrength(:,:,2));
LaS = squeeze(CcgStrength(:,:,3));
UaStrength = squeeze(CcgStrength(:,:,4));
for k = 1:size(CcgStrength,1)
    tuas = UaStrength(k,:);
    tu = UaS(k,:);
    tl = LaS(k,:);
    tinx = tu<0.0001&tl<0.0001;
    tuas(tinx) = NaN;
    UaStrength(k,:) = tuas / nansum(tuas);
end

UaStrength1 = UaStrength(som_inx,:);
figure;
plot(nansum(UaStrength1)/nansum(UaStrength1(:)),'Color','b')

UaStrength1 = UaStrength(pv_inx,:);
hold on
plot(nansum(UaStrength1)/nansum(UaStrength1(:)),'Color','r')

UaStrength1 = UaStrength(nt_inx,:);
plot(nansum(UaStrength1)/nansum(UaStrength1(:)),'Color',[0.7 0.7 0.7])