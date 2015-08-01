function aprintdata(inpdir)
%APRINTDATA   Prints figures.
%
%   See also PRINTDATA.

% Input argumnet check
error(nargchk(1,1,nargin))

% Get file list
mm = pwd;
cd(inpdir)
files = dir(inpdir);

% Open figures
for o = 3:length(files)
    fn = files(o).name;
    if ~isequal(fn(end-2:end),'fig')
        continue
    end
    open(fn)
    
% Print
    set(gcf,'PaperOrientation','landscape','PaperPositionMode','manual',...
        'PaperPosition',[0 0 28 21])
    str = ['print -f' num2str(gcf)];
    eval(str)
    close all
end
cd(mm)