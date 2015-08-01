%% F o F^(-1) = id

xx=0:0.01:8*pi;
yy=sin(xx);
[y,w] = b_fft3(yy,100);
yy2=ifft(y);
isequal(yy,yy2')
max(real(yy)-real(yy2'))
max(imag(yy)-imag(yy2'))

%% shifting the spectrum

xx=0:0.01:8*pi;
yy=sin(xx);
[y,w] = b_fft3(yy,100);
figure;plot(xx,yy)

ng = angle(y);
ng=ng-180/180*pi;
cf = abs(y) .* exp(i*ng);
cf(hlf+2:end)=flipud(conj(cf(2:hlf)));
ceeg = ifft(cf);
hold on;plot(xx,real(ceeg),'r')

%% making the spectrum "symmetric"

xx=0:0.01:8*pi;
yy=sin(xx);
[y,w] = b_fft3(yy,100);
figure;plot(xx,yy)

hlf = length(y) / 2;

y2 = y;
y2(hlf+2:end)=flipud(conj(y(2:hlf)));

iy2 = ifft(y2);
figure;plot(xx,iy2)

%% making the spectrum "symmetric"

xx=0:0.5:pi;
yy=sin(xx);
[y,w] = b_fft3(yy,100);
% figure;plot(xx,yy)

hlf = length(y) / 2;

y2 = y;
y2(hlf+1:end)=flipud(conj(y(1:hlf)));
figure;plot(real(y2));hold on;plot(real(y),'r')

%% shifting a part of the spectrum

xx=0:0.01:8*pi;
yy=sin(xx);
[y,w] = b_fft3(yy,100);
figure;plot(xx,yy)

ng = angle(y);
ng(1:500)=ng(1:500)+45/180*pi;
cf = abs(y) .* exp(i*ng);
cf(hlf+2:end)=flipud(conj(cf(2:hlf)));
ceeg = ifft(cf);
hold on;plot(xx,real(ceeg),'r')

%% two sine waves

xx=0:0.01:8*pi;
yy1=sin(xx);
yy2 = sin(xx*5);
yy=yy1+yy2;
[y,w] = b_fft3(yy,100);

figure;plot(xx,yy)
ng = angle(y);
ng(1:10)=ng(1:10)+180/180*pi;
cf = abs(y) .* exp(i*ng);
cf(hlf+2:end)=flipud(conj(cf(2:hlf)));
ceeg = ifft(cf);
hold on;plot(xx,real(ceeg),'r')

figure;plot(xx,yy)
ng = angle(y);
ng(15:25)=ng(15:25)+180/180*pi;
cf = abs(y) .* exp(i*ng);
cf(hlf+2:end)=flipud(conj(cf(2:hlf)));
ceeg = ifft(cf);
hold on;plot(xx,real(ceeg),'r')

figure;plot(xx,yy)
ng = angle(y);
ng(1:25)=ng(1:25)+180/180*pi;
cf = abs(y) .* exp(i*ng);
cf(hlf+2:end)=flipud(conj(cf(2:hlf)));
ceeg = ifft(cf);
hold on;plot(xx,real(ceeg),'r')

figure;plot(xx,yy)
ng = angle(y);
ng(1:10)=ng(1:10)+180/180*pi;
ng(15:25)=ng(15:25)+90/180*pi;
cf = abs(y) .* exp(i*ng);
cf(hlf+2:end)=flipud(conj(cf(2:hlf)));
ceeg = ifft(cf);
hold on;plot(xx,real(ceeg),'r')

%% length of data is even or odd

xx=0:0.01:8*pi;
% xx=xx(1:end-1);
yy1=sin(xx);
yy2 = sin(xx*5);
yy=yy1+yy2;
[y,w] = b_fft3(yy,100);

figure;plot(xx,yy)
ng = angle(y);
ng(1:10)=ng(1:10)+180/180*pi;
cf = abs(y) .* exp(i*ng);

leny = length(y);
mod(leny,2)
if isequal(mod(leny,2),0)
    hlf = leny / 2;
    cf(hlf+2:end)=flipud(conj(cf(2:hlf)));
else
    hlf = (leny + 1) / 2;
    cf(hlf+1:end)=flipud(conj(cf(2:hlf)));
end

ceeg = ifft(cf);
hold on;plot(xx,real(ceeg),'r')