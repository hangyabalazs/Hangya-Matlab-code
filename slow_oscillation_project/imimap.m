function gr = imimap(rI)

% rI = get(findobj(allchild(gca),'type','image'),'CData');
gr = zeros(7,9);
gr = zero2nan(gr);
for y = 2:2:8
    for x = 1:2:7
        rw = (x - 1) / 2 + 1;
        inx = (rw - 1) * 5 + y / 2;
        gr(x,y) = rI(inx,inx+1);
    end
end
for y = 1:2:9
    for x = 2:2:6
        cl = (y - 1) / 2 + 1;
        inx = cl + (x / 2 - 1) * 5;
        gr(x,y) = rI(inx,inx+5);
    end
end
% figure
% imagesc(gr,[1 2.1])