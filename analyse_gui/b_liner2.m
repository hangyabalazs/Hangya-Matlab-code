function b_liner2(plt,vs,pos,la,len,hz)
%LINER2   Plotting function for ANALYSE GUI 'Open' button.
%   LINER2(PLT,VS,POS,LA,LEN) needs five input arguments: an existing
%   line handle (PLT), a supressed data file (VS - see CREATEVS for
%   details), GUI position (POS), base of logarithm used by segmenting
%   the data in CREATEVS (LA) and length of data (LEN). It plots data
%   on ANALYSE GUI axes.
%
%   LINER2 visualizes bursts in red color. Burst data is loaded from a
%   directory specified at the beginning of the code of ANALYSE GUI.
%
%   See also ANALYSE and CREATEVS.

% Directory for burst data
global DATAPATH
global BURSTDIR
burstdir = BURSTDIR;

% Redraw data
xl = xlim;
yl = ylim;
datalim = [ceil(xl(1)*hz),ceil(xl(2)*hz)];
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
set(plt,'XData',x,'YData',vv,'Color','b');
xlim(xl);
ylim(yl);

% Visualize bursts
delete(findobj(allchild(gca),'Color','red'))

global FILENAME
fn = FILENAME;
ff = fullfile(burstdir,[fn(1:end-4) '_CLUST2.mat']);
if exist(ff,'file')
    load(ff)
else
    return
end
VBurst = vdisc(Burst);
if isequal(size(VBurst),[1,2])
    VBurst = VBurst';
end
vburst = VBurst .* (VBurst>datalim(1)) .* (VBurst<datalim(2));  % restrict vburst to the window
if isempty(vburst)
    return
end

vb1 = [vburst(1,:)>0 vburst(2,:)>0] .* [[1:size(vburst,2)] (-1)*[1:size(vburst,2)]];
if sum(vb1) > 0     % handle the case when burst exceeds the window
    ind = find(vburst,1,'last')+1;
    vburst(ind) = datalim(2);
    vb1 = [vburst(1,:)>0 vburst(2,:)>0] .* [[1:size(vburst,2)] (-1)*[1:size(vburst,2)]];
end
if sum(vb1) < 0
    ind = find(vburst,1,'first')-1;
    vburst(ind) = datalim(1);
end
fV1 = find(VBurst>datalim(2),1,'first');
fV2 = find(VBurst<datalim(1),1,'last');
if fV1 - fV2 == 1 & mod(fV1,2) == 0
    vburst = datalim;   % when one burst exceeds the window on both sides
end

vburst = vburst(find(vburst));
vburst = reshape(vburst,2,length(vburst)/2);
for k = 1:size(vburst,2)    % plot
   inx1 = length(vv) * (vburst(1,k) - datalim(1)) / (datalim(2) - datalim(1));
   inx1 = max(floor(inx1),1);
   inx2 = length(vv) * (vburst(2,k) - datalim(1)) / (datalim(2) - datalim(1));
   inx2 = min(ceil(inx2),length(vv));
   hold on
   plot(x(inx1:inx2),vv(inx1:inx2),'r')
end
hold off