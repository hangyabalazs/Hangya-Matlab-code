function nbaccuracycells_summary2
%NBACCURACYCELLS_SUMMARY2   Summary analysis of accuracy-prediting cells.
%   NBACCURACYCELLS_SUMMARY2 makes summary figures for cells predicting
%   accuracy vs. cholinergic neurons.
%
%   See also COND_ACCURACY_FR3_TEST and NBACCURACYCELLS_TEST.

% Directories
global DATAPATH
% inpdir = fullfile(DATAPATH,'NB','accuracycells','test_all',filesep);   % source files
% resdir_all = fullfile(DATAPATH,'NB','accuracycells','sign',filesep);
% resdir_ChAT = fullfile(DATAPATH,'NB','accuracycells','ChAT',filesep);
% resdir = fullfile(DATAPATH,'NB','accuracycells','summary',filesep);   % results directory
inpdir = fullfile(DATAPATH,'NB','accuracycells_newdata','test_all',filesep);   % source files
resdir_all = fullfile(DATAPATH,'NB','accuracycells_newdata','sign',filesep);
resdir_ChAT = fullfile(DATAPATH,'NB','accuracycells_newdata','ChAT',filesep);
resdir = fullfile(DATAPATH,'NB','accuracycells_newdata','summary',filesep);   % results directory
issave = true;

% Selected IDs
% nmdir = fullfile(DATAPATH,'NB','accuracycells','test_all',filesep);
nmdir = fullfile(DATAPATH,'NB','accuracycells_newdata','test_all',filesep);
files = dir(nmdir);
files2 = {files(3:end).name};
files3 = cell(size(files2));
endchar = cell2mat(regexp(files2,'_[ACFR]'));
for k = 1:length(files2)
    files3{k} = files2{k}(1:endchar(k)-1);
    inx = regexp(files3{k},'_');
    files3{k}(inx(end)) = '.';
end
cellids = unique(files3);
NumCells = length(cellids);

% p-values
[p D] = deal(nan(NumCells,6));
for k = 1:NumCells
    fnm = fullfile(inpdir,[regexprep(cellids{k},'\.','_') '_CondPerf.mat']);
    load(fnm)
    p(k,1:length(CondPerf.pValues)) = CondPerf.pValues;
    D(k,1:length(CondPerf.pValues)) = CondPerf.condDiscriminationS - ...
        CondPerf.condDiscriminationL;
end

% Significantly predictive cells
inx = any(p<0.01&D>0,2);
main(cellids(inx),inpdir,resdir_all,issave)

% Cholinergic cells
ChAT = selectcell(['"ChAT+"==1&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("session_type",{''gonogo'',''feedbackdelay''})&' ...
    'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})']);  % identified ChAT+ cells
ChAT = [ChAT 'n045_121217x_4.6'];   % clustered based on light-evoked spikes
pChAT = selectcell(['"pChAT+"==1&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})']);  % putative ChAT+ cells
% allChAT = [ChAT pChAT];
allChAT = ChAT;
main(ChAT,inpdir,resdir_ChAT,issave)
% main(allChAT,inpdir)

% -------------------------------------------------------------------------
function main(cellids,inpdir,resdir,issave)

% Difference of discrimination in high - low FR trials
NumCells = length(cellids);   % number of NB cells
Discrimination1 = nan(NumCells,6);
Discrimination2 = nan(NumCells,3);
for iC = 1:NumCells   % loop through NB cells
    cellid = cellids{iC};
    cellidt = regexprep(cellid,'\.','_');
    fnm = [inpdir cellidt '_CondPerf.mat'];
    load(fnm)
    df = CondPerf.condDiscriminationS - CondPerf.condDiscriminationL;
    df = df(~isnan(df));
    if isempty(df)
        Discrimination2(iC,1:3) = [NaN NaN NaN];
    elseif length(df) == 1
        Discrimination2(iC,1:3) = [NaN NaN df];
    else
        Discrimination2(iC,1) = df(1);
        Discrimination2(iC,2) = nanmean(df(2:end-1));
        Discrimination2(iC,3) = df(end);
    end
    Discrimination1(iC,1:length(CondPerf.pValues)) = ...
        CondPerf.condDiscriminationS - CondPerf.condDiscriminationL;  % difference of discrimination in high - low FR trials
    if issave
        fnm = [inpdir cellidt '_FRDPrime.jpg'];
        copyfile(fnm,resdir)
    end
end

% Plot
Discrimination1(all(isnan(Discrimination1),2),:) = [];
figure
imagesc(Discrimination1)

figure
plot(nanmean(Discrimination1))

figure
% plot(nanmean(Discrimination2))
errorshade(1:3,nanmean(Discrimination2),nanse(Discrimination2),...
    'LineColor','k','ShadeColor','k')
keyboard