function vm = b_rvm(n,mu,k)
%RVM   Generates random numbers from a von Mises distribution.
%   VM = RVM(N,MU,K) generates N length independent random sample from a
%   von Mises distribution with MU mean and K concentratzion parameter.
%
%   See also RMIXEDVM and RRMIXEDVM.

% Input argument check
error(nargchk(3,3,nargin))

% Random von Mises sample
vm = (1:n);
a = 1 + (1 + 4 * (k^2))^0.5;
b = (a - (2 * a)^0.5) / (2 * k);
r = (1 + b^2) / (2 * b);
obs = 1;
while obs <= n
    U1 = rand(1);
    z = cos(pi*U1);
    f = (1 + r * z) / (r + z);
    c = k * (r - f);
    U2 = rand(1);
    if c * (2 - c) - U2 > 0
        U3 = rand(1);
        vm(obs) = sign(U3-0.5) * acos(f) + mu;
        vm(obs) = mod(vm(obs),2*pi);
        obs = obs + 1;
    elseif log(c/U2) + 1 - c >= 0
            U3 = rand(1);
            vm(obs) = sign(U3-0.5) * acos(f) + mu;
            vm(obs) = mod(vm(obs),2*pi);
            obs = obs + 1;
    end
end