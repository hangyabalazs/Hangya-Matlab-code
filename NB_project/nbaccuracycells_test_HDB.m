function nbaccuracycells_test_HDB
%NBACCURACYCELLS_TEST   Cells predicting accuracy.
%   NBACCURACYCELLS_TEST performs analysis to asess whether cells are
%   related to accuracy. See details about the analysis in
%   COND_ACCURACY_FR3_TEST.
%
%   See also COND_ACCURACY_FR3_TEST.

%   Balazs Hangya, Cold Spring Harbor Laboratory
%   1 Bungtown Road, Cold Spring Harbor
%   balazs.cshl@gmail.com
%   3-Dec-2013

%   Edit log: BH 12/3/13

% Directories
global DATAPATH
% resdir = fullfile(DATAPATH,'NB','accuracycells','test_good',filesep);   % results directory
resdir = fullfile(DATAPATH,'HDB','accuracycells_newdata','test_all',filesep);   % results directory
issave = true;

% Cells
% nmdir = fullfile(DATAPATH,'NB','accuracycells','good');
nmdir = fullfile(DATAPATH,'HDB','accuracycells_newdata','view');
files = dir(nmdir);
files2 = {files(3:end).name};
files3 = cell(size(files2));
endchar = cell2mat(regexp(files2,'_F'));
for k = 1:length(files2)
    files3{k} = files2{k}(1:endchar(k)-1);
    inx = regexp(files3{k},'_');
    files3{k}(inx(end)) = '.';
end
cellids = unique(files3);
NumCells = length(cellids);   % number of NB cells
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