function b_liner2_sniffwhisk(plt,vs,pos,la,len,hz,fp)
%LINER2_SNIFFWHISK   Plotting function for ANALYSE GUI 'Open' button.
%   LINER2_SNIFFWHISK(PLT,VS,POS,LA,LEN,FP) needs six input arguments: an
%   existing line handle (PLT), a supressed data file (VS - see CREATEVS
%   for details), GUI position (POS), base of logarithm used by segmenting
%   the data in CREATEVS (LA), length of data (LEN) and first time stamp
%   (FP). It plots data on ANALYSE GUI axes.
%
%   LINER2_SNIFFWHISK visualizes bursts in red color. Burst data is loaded
%   from a directory specified at the beginning of the code of ANALYSE GUI.
%
%   LINER2_SNIFFWHISK plots Neuralynx events.
%
%   See also ANALYSE_SNIFFWHISK and CREATEVS.

% Directory for burst data
% global DATAPATH
% global BURSTDIR
% burstdir = BURSTDIR;

% Redraw data
xl = get(get(get(0,'CurrentFigure'),'CurrentAxes'),'XLim');
yl = get(get(get(0,'CurrentFigure'),'CurrentAxes'),'YLim');
datalim = [ceil(xl(1)*hz),ceil(xl(2)*hz)] - fp * hz;
width = pos(3);
numdata = width / (datalim(2) - datalim(1) + eps) * len * 2;

if numdata >= len
    vvind = 1;
else
    vvind = max(1,length(vs)-ceil(log(numdata)/log(la)));
end
vv = vs{vvind};
ind1 = max(round(datalim(1)*length(vv)/len)+1,1);
ind2 = min(round(datalim(2)*length(vv)/len)-1,length(vv));
vv = vv(ind1:ind2);
x = [];
if ~isempty(vv)
    x = linspace(xl(1), xl(2), length(vv));
end
set(plt,'XData',x,'YData',vv,'Color','k');
set(get(get(0,'CurrentFigure'),'CurrentAxes'),'XLim',xl);
set(get(get(0,'CurrentFigure'),'CurrentAxes'),'YLim',yl);