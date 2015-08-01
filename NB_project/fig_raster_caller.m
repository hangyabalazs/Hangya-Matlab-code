% ChAT cells

% selstr = ['"ChAT+"==1&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
%     'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})'];
% selstr = '"ChAT+"==1&"ID_PC">20&"Lr_PC"<0.15&"validity"==1';
% ChAT = selectcell(selstr);   % cell IDs for ChAT cells
% ChAT = [ChAT 'n067_141017a_1.3'];   % without miss-triggering it passes the cluster criteria
% ChAT = [ChAT 'n067_141019a_5.2'];   % light spike assisted clustering

ChAT = selectcell(['"ChAT+"==1&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("session_type",{''gonogo'',''feedbackdelay''})&' ...
    'ismember("Area1",{''HDB'',''SI'',''VP''})']);  % identified (n = 12)
ChAT = [ChAT 'n067_141017a_1.3'];   % without miss-triggering it passes the cluster criteria
ChAT = [ChAT 'n067_141019a_5.2'];   % light spike assisted clustering

% Directories
global DATAPATH
fs = filesep;
resdir = [DATAPATH 'HDB' fs 'raster_PSTH' fs];

NumChAT = length(ChAT);
for iC = 1:NumChAT
    cellid = ChAT{iC};
    fig_raster(cellid)
    H = gcf;
    fnm = fullfile(resdir,['RASTER_' regexprep(cellid,'\.','_') '.fig']);
    saveas(H,fnm)
    
    % Print
    set(H,'PaperPositionMode','auto')
    orient tall
    str = ['print -f' num2str(H)];
    eval(str)
    close(H)
end