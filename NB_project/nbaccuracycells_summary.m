function nbaccuracycells_summary
%NBACCURACYCELLS_SUMMARY   Summary analysis of accuracy-prediting cells.
%   NBACCURACYCELLS_SUMMARY makes summary figures for cells predicting
%   accuracy vs. cholinergic neurons.
%
%   See also COND_ACCURACY_FR3_TEST and NBACCURACYCELLS_TEST.

% Directories
global DATAPATH
inpdir = fullfile(DATAPATH,'NB','accuracycells',filesep);   % source files
resdir = fullfile(DATAPATH,'NB','accuracycells','summary',filesep);   % results directory
issave = true;

% Selected cells
nmdir = fullfile(DATAPATH,'NB','accuracycells','good2',filesep);
files = dir(nmdir);
files = cellfun(@(s)[s(1:14) '.' s(16)],{files(3:end).name},...
    'UniformOutput',false);
cellids = unique(files);
main(cellids,inpdir)

% Cholinergic cells
ChAT = selectcell(['"ChAT+"==1&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})']);  % identified ChAT+ cells
ChAT = [ChAT 'n045_121217x_4.6'];   % clustered based on light-evoked spikes
pChAT = selectcell(['"pChAT+"==1&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})']);  % putative ChAT+ cells
allChAT = [ChAT pChAT];
main(ChAT,inpdir)
main(allChAT,inpdir)

% -------------------------------------------------------------------------
function main(cellids,inpdir)

% Difference of discrimination in high - low FR trials
NumCells = length(cellids);   % number of NB cells
Discrimination = nan(NumCells,4);
for iC = 1:NumCells   % loop through NB cells
    cellid = cellids{iC};
    cellidt = regexprep(cellid,'\.','_');
    fnm = [inpdir cellidt '_CondPerf.mat'];
    load(fnm)
    Discrimination(iC,:) = CondPerf.condDiscriminationS - CondPerf.condDiscriminationL;  % difference of discrimination in high - low FR trials
end

% Plot
Discrimination(all(isnan(Discrimination),2),:) = [];
figure
imagesc(Discrimination)

figure
plot(nanmean(Discrimination))
keyboard