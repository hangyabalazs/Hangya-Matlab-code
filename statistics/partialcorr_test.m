% Test program for partial correlation.

% SCENARIO1: Rxy_c_z ~ Rxy
% SCENARIO2: Rxy_c_z ~ 0
% SCENARIO3: Rxy_c_z < Rxy
% SCENARIO4: Rxy_c_z ~ 0

% X = rand(1000,1);   % phase
% Y = 2*X;            % rt
% Z = rand(1000,1);   % amp
% Z = X / 2;;

X = rand(1000,1);   % SCENARIO1: RT and amp are dependent on phase
Z = X + 2*rand(1000,1);
Y = X + 2*rand(1000,1);

% X = rand(1000,1);   % scenario: RT and amp are dependent on phase, RT partially depends on amp 
% Z = X/3 + 3*rand(1000,1);
% Y = X + 3*rand(1000,1) + 50*Z;

% X = rand(1000,1);   % SCENARIO2: RT is dependent on phase, amp is dependent on phase and RT 
% Y = X + 2*rand(1000,1);
% Z = X + Y + 2*rand(1000,1);

% X = rand(1000,1);   % SCENARIO3: amp is dependent on phase, RT is dependent on phase and amp 
% Z = X + 2*rand(1000,1);
% Y = X + Z + 2*rand(1000,1);

% X = rand(1000,1);   % SCENARIO4: amp is dependent on phase, RT is dependent on amp 
% Z = X + 2*rand(1000,1);
% Y = Z + 2*rand(1000,1);

X1 = [ones(size(X,1),1) X];
Y1 = [ones(size(Y,1),1) Y];
Z1 = [ones(size(Z,1),1) Z];
[b,bint,rY,rint,regstats] = regress(Y,Z1);
[b,bint,rX,rint,regstats] = regress(X,Z1);

corrcoef(rY,rX)

% -------------------------------------------------------------------------

pR = corrcoef(X,Y);
Rxy = pR(2);
pR = corrcoef(X,Z);
Rxz = pR(2);
pR = corrcoef(Y,Z);
Ryz = pR(2);

Rxy_c_z = (Rxy - Rxz * Ryz) / (sqrt(1-Rxz^2) * (sqrt(1-Ryz^2)))

Rxy
Rxy_c_z / Rxy
