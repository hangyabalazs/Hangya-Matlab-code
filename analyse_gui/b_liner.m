function b_liner(plt,vs,pos,la,len,hz)
%LINER   Plotting function for ANALYSE GUI 'Open' button.
%   LINER(PLT,VS,POS,LA,LEN,HZ) needs six input arguments: an existing
%   line handle (PLT), a supressed data file (VS - see CREATEVS for
%   details), GUI position (POS), base of logarithm used by segmenting
%   the data in CREATEVS (LA), length of data (LEN) and sampling rate (HZ).
%   It plots data on ANALYSE GUI axes.
%
%   See also ANALYSE and CREATEVS.

xl = xlim;
yl = ylim;
datalim = [ceil(xl(1)*hz),ceil(xl(2)*hz)];
width = pos(3);
numdata = width / (datalim(2) - datalim(1) + eps) * len * 2;

if numdata >= len
    vv = vs{1};
else
    vv = vs{max(1,end-ceil(log(numdata)/log(la)))};
end
ind1 = max(round(datalim(1)*length(vv)/len)+1,1);
ind2 = min(round(datalim(2)*length(vv)/len)-1,length(vv));
vv = vv(ind1:ind2);
x = [];
if ~isempty(vv)
    x = linspace(xl(1), xl(2), length(vv));
end
set(plt,'XData',x,'YData',vv,'Color','b');
xlim(xl);
ylim(yl);