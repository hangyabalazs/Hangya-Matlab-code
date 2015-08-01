%%

a = 5;
b = 2;
c = 288;
d = 2945-288;


%%

p = nchoosek(a+b,a) * nchoosek(c+d,c) / nchoosek(a+b+c+d,a+c)

%%

[h p stat] = b_chi2test2([a b],[c d])