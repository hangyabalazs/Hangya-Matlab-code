% Import

I=imread('F:\balazs\_personal\kepek\Obergurgl_2008\n640434547_255012_6671.jpg');
% figure;imagesc(I)
% I1=cat(3,I(:,:,1),zeros(453,604),zeros(453,604))
% figure;imagesc(I1)
Ired=I(:,:,1);
Ired=double(Ired);
% figure;imagesc(Ired)
red1=Ired(1,:);
red2=Ired(2,:);
red10 = Ired(10,:);
figure;plot(red1)
hold on;plot(red2,'r')
m1 = red1;
m2 = red10;

fr = fir1(64,10/300,'high');
fred300=filtfilt(fr,1,m2);
figure;plot(fred300)
vred=disc(fred300,4);

% ffred1=fft(red1);
% figure;plot(ffred1)
% figure;plot(ffred1.*conj(ffred1))
% pred1=ffred1.*conj(ffred1);
% figure;plot(pred1(2:300))
% figure;plot(pred1(10:300))

fr=fir1(64,[10/300 40/300]);
fred1=filtfilt(fr,1,m1);
figure;plot(fred1)
afre=angle(hilbert(fred1));
ang=afre(vred);
b_mvl(ang)
figure;rose(ang)