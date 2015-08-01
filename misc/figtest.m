H=figure
for i = 1:10
    figure(H)
    set(H,'HitTest','off')
    plot([i i])
    pause(5)
end