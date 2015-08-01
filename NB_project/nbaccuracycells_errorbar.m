function nbaccuracycells_errorbar
%NBACCURACYCELLS_ERRORBAR   Cells predicting accuracy.
%   NBACCURACYCELLS_ERRORBAR performs analysis to asess whether cells are
%   related to accuracy. See details about the analysis in
%   COND_ACCURACY_FR3_TEST.
%
%   See also COND_ACCURACY_FR3_ERRORBAR.

%   Balazs Hangya, Cold Spring Harbor Laboratory
%   1 Bungtown Road, Cold Spring Harbor
%   balazs.cshl@gmail.com
%   3-Dec-2013

%   Edit log: BH 12/3/13, 1/11/14

% Directories
global DATAPATH
% inpdir = fullfile(DATAPATH,'NB','accuracycells','test_all',filesep);   % source files
% resdir = fullfile(DATAPATH,'NB','accuracycells','errorbar',filesep);   % results directory
inpdir = fullfile(DATAPATH,'HDB','accuracycells_newdata','test_all',filesep);   % source files
resdir = fullfile(DATAPATH,'HDB','accuracycells_newdata','errorbar',filesep);   % results directory
issave = true;

% Cells
% Selected IDs
files = dir(inpdir);
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
cellids = cellids(inx);
NumCells = length(cellids);

for iC = 1:NumCells   % loop through NB cells
    cellid = cellids{iC};
    disp(cellid)
    
    % Skip if no behavior
    ibh = getvalue('session_type',cellid);
    if isequal(ibh{1},'no behavior')
        continue
    end
    
    % Accuracy cells, action cells
    [H CondPerf] = cond_accuracy_fr3_test(cellid,'window',[-0.5 0],'limit2iti',true,'event','StimulusOn','display',true);
    cellidt = regexprep(cellid,'\.','_');
    if issave   % save
        fnm = [resdir cellidt '_CondPerf.mat'];
        save(fnm,'CondPerf')
        
        fnm = [resdir cellidt '_FRDPrime.fig'];
        saveas(H.Hdpr,fnm)
        fnm = [resdir cellidt '_FRDPrime.jpg'];
        saveas(H.Hdpr,fnm)
    end
    close all
end