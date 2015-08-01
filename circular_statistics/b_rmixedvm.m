function R = b_rmixedvm(n,mu1,mu2,kappa1,kappa2,p)
%RMIXEDVM   Generates random numbers from a mixture of two von Mises distributions.
%   R = RMIXEDVM(N,MU1,MU2,KAPPA1,KAPPA2,P) generates N length independent
%   random mixed von Mises sample from mixture of MU1 mean KAPPA1
%   concentration and MU2 mean KAPPA2 concentration von Mises distributions
%   using P mixing ratio. Random numbers are from density  function
%   p*VM(mu1, kappa1) + (1-p)*VM(mu2, kappa2) where VM is the von Mises 
%   density function.
%
%   See also RVM and RRMIXEDVM.

% Input argument check
error(nargchk(6,6,nargin))

% Random mixed von Mises sample
R = zeros(1,n);
for i = 1:n
    test = rand(1);
    if test < p
        R(i) = b_rvm(1,mu1,kappa1);
    else
        R(i) = b_rvm(1,mu2,kappa2);
    end
end