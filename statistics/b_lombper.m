function [pxx,pyy,jmax,prob,z,effm,period_handle]=b_lombper(x,y)
%LOMBPER   Lomb-Scargle periodogram.
%   LOMBPER(X,Y) calculates Lomb-Scargle periodogram where X contains the sampling points
%   and Y contains the sampled data.
%
%   See also LOMB_PERIOD and PERIODRUN.


% best values ofac=4, hifac=1
ofac=16; hifac=4;

% disp(' calculation of Lomb - Scargle periodogram ')

% fsamp=input(' sample frequency [m] = ');
% fnyq=fsamp/2;

wb = waitbar(0,'Please wait...','Position',[360 250 275 50]);    %Progress indicator

n=length(x);

wi=zeros(size(1:n));
wpi=zeros(size(1:n));
wpr=zeros(size(1:n));
wr=zeros(size(1:n));



nout=.5*ofac*hifac*n;

py=zeros(size(1:nout));
px=zeros(size(1:nout));

avey=mean(y);
vary=cov(y);

xmax=max(x);
xmin=min(x);
xdif=xmax-xmin;
xave=.5*(xmax+xmin);

pymax=0.;

pnow=1./(xdif*ofac);

for j=1:n
    arg=2*pi*((x(j)-xave)*pnow);
    wpr(j)=(-2.)*(sin(.5*arg))^2;
    wpi(j)=sin(arg);
    wr(j)=cos(arg);
    wi(j)=wpi(j);
end

for i=1:nout
    px(i)=pnow;
    sumsh=0.;
    sumc=0.;
    
    for j=1:n
        c=wr(j);
        s=wi(j);
        sumsh=sumsh+s*c;
        sumc=sumc+(c-s)*(c+s);
    end
    
    wtau=.5*atan2(2.*sumsh,sumc);
    swtau=sin(wtau);
    cwtau=cos(wtau);
    sums=0.;
    sumc=0.;
    sumsy=0.;
    sumcy=0.;
    for j=1:n
        s=wi(j);
        c=wr(j);
        ss=s*cwtau-c*swtau;
        cc=c*cwtau+s*swtau;
        sums=sums+ss*ss;
        sumc=sumc+cc*cc;
        yy=y(j)-avey;
        sumsy=sumsy+yy*ss;
        sumcy=sumcy+yy*cc;
        wtemp=wr(j);
        wr(j)=(wr(j)*wpr(j)-wi(j)*wpi(j))+wr(j);
        wi(j)=(wi(j)*wpr(j)+wtemp*wpi(j))+wi(j);
        % bno( 'at the end of the last inner loop '),keyboard
    end
    % bno( 'before calc. py(i) '),keyboard
    
    py(i)=.5*(sumcy*sumcy/sumc+sumsy*sumsy/sums)/vary;
    if py(i)>=pymax
        jmax=i;
        pymax=py(i);
    end
    pnow=pnow+1./(ofac*xdif);
    if mod(i,100)==0
        waitbar(i/nout)   %Progress indicator
    end
end

close(wb);   %Close progress indicator

expy=exp(-pymax);
effm=2.*nout/ofac;
prob=effm*expy;
if(prob>0.01)
    prob=1.0-(1.0-expy)^effm;
end


clear wr wpr wpi wi

fnyq=n/(2*xdif);
% disp([ ' Nyquist frequency = ' num2str(fnyq)])

% significance levels
% for p=.05, .02, .01 .001

z05=log(effm/.05);
z02=log(effm/.02);
z01=log(effm/.01);
% z31=log(effm/.1586555);

% for p=.001
z001=log(1)-log(1-exp((log(1-.001))/effm));

% z=[ z31 ];
z=[ z05 z02 z01 z001 ];


% plot Lomb-spectrum and signif. levels

% px1=px*1000;
px1=px;
llpx=length(px1);
pxx=[0 px1 px1(llpx) ];
pyy=[0 py 0];

period_handle = fill(pxx,pyy,'r');
y_lim = ylim;
axis([0 max(pxx) y_lim(1) y_lim(2)])
hold on
v=axis;
% plot([v(1) v(2)], [ z31 z31], 'k:')
% text(v(2),z31,' 31% ')
plot([v(1) v(2)], [ z05 z05], 'k:')
text(v(2),z05,' 5% ')
plot([v(1) v(2)], [ z02 z02], 'k:')
text(v(2)+v(2)*.03,z02,' 2% ')
plot([v(1) v(2)], [ z01 z01], 'k:')
text(v(2)+v(2)*.06,z01,' 1% ')
plot([v(1) v(2)], [ z001 z001], 'k:')
text(v(2),z001,' .1% ')
xlabel(' frequency [Hz]')
ylabel(' PSD [s^2]')
% title(' Lomb - Scargle periodogram ')
hold off