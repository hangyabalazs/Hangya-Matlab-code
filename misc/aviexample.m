%% AVI example

fig=figure;
set(fig,'DoubleBuffer','on');
set(gca,'xlim',[-80 80],'ylim',[-80 80],...
    'NextPlot','replace','Visible','off')
mov = avifile('example.avi')
x = -pi:.1:pi;
radius = 0:length(x);
for k=1:length(x)
 h = patch(sin(x)*radius(k),cos(x)*radius(k),...
    [abs(cos(x(k))) 0 0]);
 set(h,'EraseMode','xor');
 F = getframe(gca);
 mov = addframe(mov,F);
end
mov = close(mov);

%% Movie


for k = 1:16
 plot(fft(eye(k+16)))
 axis equal
 M(k) = getframe;
end

movie(M,30)
movie2avi(M,'example3.avi','Compression','none')

%% Mpeg

c = colormap;
mpgwrite(M,c,'proba')