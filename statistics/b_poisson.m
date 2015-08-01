function b_poisson(lambda)
%POISSON   Poisson density function.
%   POISSON(LAMBDA) plots Poisson density function with LAMBDA parameter.

xx = [1:100];
for i = 1:100
    yy(i) = (lambda ^ xx(i)) / (factorial(xx(i))) * exp(-lambda);
end
plot(xx,yy)