%%

yd=get(gco,'YData');
ydd=yd;
xd=get(gco,'XData');
inx=yd==0;
ydd(inx)=[];
xdd=xd;
xdd(inx)=[];
% figure;plot(xdd,ydd)
yds = ydd;
xds = xdd;

%%

yd=get(gco,'YData');
ydd=yd;
xd=get(gco,'XData');
inx=yd==0;
ydd(inx)=[];
xdd=xd;
xdd(inx)=[];
yds = [yds; ydd];
xds = [xds; xdd];

%%

arls = arrowlength;
alxp = alx;

%%

isequal(alxp,alx)

%%

arls = [arls; arrowlength];