function b_animation
%ANIMATION   Demo program for animation.
%   This program shows how to program animation in Matlab through 
%   demonstrating Brownian motion. You are able to change the number
%   of moving points, the velocity and the pause time in the first
%   three code lines.

n = 20;    %number of points
s = 0.006;  %velocity
k = 0.1;   %pause time
disp(['To quit animation, press Ctrl+c!'])
x=rand(n,1)-0.5;
y=rand(n,1)-0.5;
h = plot(x,y,'.');
title(['Brownian motion'])
axis([-1 1 -1 1 ]);
axis square
grid off
set(h,'EraseMode','xor','MarkerSize',18)
try
    while 1
        drawnow
        x = x+s*randn(n,1);
        y = y+s*randn(n,1);
        set(h,'XData',x,'YData',y)
        pause(k)
    end
catch
    return
end