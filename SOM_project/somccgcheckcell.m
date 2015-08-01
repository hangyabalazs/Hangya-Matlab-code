%%

cnt = 0;
for k = 1:length(nc1)
    if ~any(nc2>nc1(k)+20&nc2<nc1(k)+20+5)
        cnt = cnt + 1;
    end
end

%%

cnt = zeros(1,20);
for t = 1:20
    cnt(t) = 0;
    for k = 1:length(nc1)
        cnt(t) = cnt(t) + length(find(nc2>nc1(k)+t-1&nc2<=nc1(k)+t-1+3));
    end
end
figure
plot(cnt)

%%

cnt = zeros(1,20);
for t = 1:20
    cnt(t) = 0;
    for k = 1:length(nc2)
        cnt(t) = cnt(t) + length(find(nc1>nc2(k)+t-1&nc1<=nc2(k)+t-1+3));
    end
end
figure
plot(cnt)