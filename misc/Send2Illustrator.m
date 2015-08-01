% Function to open MATLAB plot in Adobe Illustrator
% JS 2013
function Send2Illustrator(FigureHandle)
UP = userpath;
UP = [UP(1:(end-1)) '\'];
TempFileDir = [UP 'IllustratorTemp\'];
if ~isdir(TempFileDir)
    mkdir(TempFileDir);
end
TempFilePath = [TempFileDir 'IllustratorTemp.eps'];
IllustratorPath = '"C:\Program Files (x86)\Adobe\Adobe Illustrator CS5\Support Files\Contents\Windows\Illustrator.exe"';
OriginalRenderer = get(FigureHandle, 'Renderer');
set(FigureHandle, 'Renderer', 'painters');
% Set all fonts to arial 12pt
ChildPlots = get(FigureHandle, 'Children');
nPlots = length(ChildPlots);
for x = 1:nPlots
    CurrentPlot = ChildPlots(x);
    PlotType = get(CurrentPlot, 'Type');
    if strcmp(PlotType, 'axes')
        set(CurrentPlot, 'FontName', 'Arial');
        CurrentFontSize = get(CurrentPlot, 'FontSize');
        if CurrentFontSize < 12
            set(CurrentPlot, 'FontSize', 12);
        end
        XLabelHandle = get(CurrentPlot, 'XLabel');
        set(XLabelHandle, 'FontName', 'Arial');
        CurrentFontSize = get(XLabelHandle, 'FontSize');
        if CurrentFontSize < 12
            set(XLabelHandle, 'FontSize', 12);
        end
        YLabelHandle = get(CurrentPlot, 'YLabel');
        set(YLabelHandle, 'FontName', 'Arial');
        CurrentFontSize = get(YLabelHandle, 'FontSize');
        if CurrentFontSize < 12
            set(YLabelHandle, 'FontSize', 12);
        end
    end
end
set(FigureHandle,'PaperPositionMode','auto')
print(FigureHandle, TempFilePath, '-depsc');


fid = fopen(TempFilePath);
ff = char(fread(fid))';
fclose(fid);
figure(FigureHandle);
actualfont = get(gca,'FontName');
mlabfontlist = {'AvantGarde','Helvetica-Narrow','Times-Roman','Bookman',...
    'NewCenturySchlbk','ZapfChancery','Courier','Palatino','ZapfDingbats',...
    'Helvetica'};
for k = 1:length(mlabfontlist)
ff = strrep(ff,mlabfontlist{k},actualfont);
end
% open the file up and overwrite it
fid = fopen(TempFilePath,'w');
fprintf(fid,'%s',ff);
fclose(fid);

set(FigureHandle, 'Renderer', OriginalRenderer);
if ispc
    IllustratorLaunchCommand = ['Start ' IllustratorPath ' ' TempFilePath];
else
    IllustratorLaunchCommand = [IllustratorPath ' ' TempFilePath];
end
system(IllustratorLaunchCommand);