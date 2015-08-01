% ChAT cells

ChAT = {'n071_141129a_2.3' 'n071_141201a_2.3' 'n071_141201a_7.2' 'n071_141218a_6.1'...
    'n071_141229a_3.2' 'n072_141222a_4.3' 'n072_141223a_4.2' 'n072_150103a_2.4'};

% Directories
global DATAPATH
fs = filesep;
resdir = [DATAPATH 'NB' fs 'raster_PSTH_newdata' fs];

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