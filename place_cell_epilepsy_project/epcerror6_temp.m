function fc = epcerror6(xc)
%EPCERROR6   Error function for EPCFIT4.
%   FC = EPCERROR6(XC) returns least square error in FC. Input argument XC
%   should be 1-by-3 array containing constant 'B', constant 'A' and
%   constant 'M' parameters for EPCMETAMODEL6.
%
%   FMINSEARCHBND calls EPCERROR6 for minimization.
%
%   See also FMINSEARCHBND, EPCMETAMODEL6 and EPCFIT6.

% Get parameters
constB = xc(1);
constA = xc(2);
constM = xc(3);
k_LP0 = xc(4);

% Call the model
[PC LC T] = epcmetamodel6(constB,constA,constM,k_LP0);

% Count the iterations
global COUNTER
COUNTER = COUNTER + 1;

% Get the series of check points
cps = 180 + [-120 -60 30 3*60 3*24*60-120 3*24*60-60 3*24*60+30 ...
    3*24*60+3*60 7*24*60-120 7*24*60-60 7*24*60+30 7*24*60+3*60];
pcs = arrayfun(@(s)PC(find(T<s,1,'last')),cps);
lcs = arrayfun(@(s)LC(find(T<s,1,'last')),cps);

% Real data
realPC = [77.36 70.75 53.85 59.80 52.17 51.79 39.81 41.67 31.63 27.66 13.48 30.12];

% Least square error
R = sum((pcs-realPC).^2);
disp([COUNTER R])

% Output
fc = R;
global DYNAMICS
DYNAMICS.CONSTB = constB;
DYNAMICS.CONSTA = constA;
DYNAMICS.CONSTM = constM;
DYNAMICS.K_LP0 = k_LP0;
DYNAMICS.PC = PC;
DYNAMICS.LC = LC;
DYNAMICS.T = T;
DYNAMICS.REALTIME = cps;
DYNAMICS.REALDATA = realPC;
DYNAMICS.PCS = pcs;
DYNAMICS.LCS = lcs;
DYNAMICS.ERROR = R;