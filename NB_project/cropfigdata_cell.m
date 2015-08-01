%% crop figures beyond the axes before export to eps

xd=get(gco,'xdata');
yd=get(gco,'ydata');
% inx=xd<-.6|xd>.6;
inx=xd>0;
xd(inx)=[];
yd(inx)=[];
set(gco,'xdata',xd,'ydata',yd)

%% images

tmp=get(gco,'cdata');
xd=get(gco,'xdata');
% inx=xd<-.6|xd>.6;
inx=xd>0;
xd=get(gco,'xdata');
tmp(:,inx)=[];
xd(inx)=[];
set(gco,'cdata',tmp,'xdata',xd)