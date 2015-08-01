function nbaccuracycells_HDB
%NBACCURACYCELLS   Cells predicting accuracy.
%   NBACCURACYCELLS performs analysis to asess whether cells are related to
%   accuracy. See details about the analysis in COND_ACCURACY_FR3.
%
%   See also COND_ACCURACY_FR3.

%   Balazs Hangya, Cold Spring Harbor Laboratory
%   1 Bungtown Road, Cold Spring Harbor
%   balazs.cshl@gmail.com
%   20-Nov-2013

%   Edit log: BH 11/20/13

% Directories
global DATAPATH
resdir = fullfile(DATAPATH,'HDB','accuracycells_newdata',filesep);   % results directory
issave = true;

% Areas
NB = selectcell(['"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("session_type",{''gonogo'',''feedbackdelay''})&' ...
    'ismember("Area1",{''HDB'',''SI'',''VP''})']);

% Cell types
ChAT = selectcell(['"ChAT+"==1&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("session_type",{''gonogo'',''feedbackdelay''})&' ...
    'ismember("Area1",{''HDB'',''SI'',''VP''})']);  % identified ChAT+ cells
ChAT = [ChAT 'n067_141017a_1.3'];   % without miss-triggering it passes the cluster criteria
ChAT = [ChAT 'n067_141019a_5.2'];   % light spike assisted clustering
pChAT = selectcell(['"pChAT+"==1&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})']);  % putative ChAT+ cells
PV = selectcell(['"PV+"==1&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})']);  % identified PV+ cells

NumCells = length(NB);   % number of NB cells
for iC = 1:NumCells   % loop through NB cells
    cellid = NB{iC};
    
    % Skip if no behavior
    ibh = getvalue('session_type',cellid);
    if isequal(ibh{1},'no behavior')
        continue
    end
    
    % Accuracy cells, action cells
    [H CondPerf] = cond_accuracy_fr3(cellid,'window',[-0.5 0],'limit2iti',true,'event','StimulusOn','display',true);
    cellidt = regexprep(cellid,'\.','_');
    if issave   % save
        fnm = [resdir cellidt '_CondPerf.mat'];
        save(fnm,'CondPerf')
        
        fnm = [resdir cellidt '_AccuracyCRT.fig'];
        saveas(H.Hca,fnm)
        fnm = [resdir cellidt '_AccuracyCRT.jpg'];
        saveas(H.Hca,fnm)
        
        fnm = [resdir cellidt '_ResponseCRT.fig'];
        saveas(H.Hcr,fnm)
        fnm = [resdir cellidt '_ResponseCRT.jpg'];
        saveas(H.Hcr,fnm)
        
        fnm = [resdir cellidt '_RTAccuracy.fig'];
        saveas(H.Hrta,fnm)
        fnm = [resdir cellidt '_RTAccuracy.jpg'];
        saveas(H.Hrta,fnm)
        
        fnm = [resdir cellidt '_FRAccuracy.fig'];
        saveas(H.Hfra,fnm)
        fnm = [resdir cellidt '_FRAccuracy.jpg'];
        saveas(H.Hfra,fnm)
        
        fnm = [resdir cellidt '_FRResponse.fig'];
        saveas(H.Hfrr,fnm)
        fnm = [resdir cellidt '_FRResponse.jpg'];
        saveas(H.Hfrr,fnm)
        
        fnm = [resdir cellidt '_FRDPrime.fig'];
        saveas(H.Hdpr,fnm)
        fnm = [resdir cellidt '_FRDPrime.jpg'];
        saveas(H.Hdpr,fnm)
        
        fnm = [resdir cellidt '_CondResp.mat'];
        save(fnm,'CondPerf')
    end
    close all
end