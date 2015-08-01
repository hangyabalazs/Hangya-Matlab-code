%%

tic
mu = 0.1;
r = -mu .* log(rand(1,1000));
pp = cumsum(r);
pp(pp>50) = [];
toc

%%

tic
mu = 0.1;
pp = 0;
while pp(end) < 50
    pp = [pp pp(end)-mu .* log(rand)];
end
toc