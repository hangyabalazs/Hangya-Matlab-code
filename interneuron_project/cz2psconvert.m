%%

ux = u2 / max(u2) * mx;
bs = [280:20:340 0:20:260];

ux2 = ux(1:18)';
cuvs = exp(i*(bs/180*pi));
wcuvs = cuvs .* ux2;
mv = sum(wcuvs);
mang = angle(mv);
mangdeg = mang * 180 / pi;
