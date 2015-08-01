function fc = epcerror(xc)
%EPCERROR   Error function for EPCFIT.
%   FC = EPCERROR(XC) returns least square error in FC. Input argument XC
%   should be 1-by-2 array containing constant 'B' and constant 'A'
%   parameters for EPCMETAMODEL.
%
%   FMINSEARCHBND calls EPCERROR for minimization.
%
%   See also FMINSEARCHBND, EPCMETAMODEL and EPCFIT.

% Get parameters
constB = xc(1);
constA = xc(2);
k_PL0 = xc(3);

% Call the model
[PC LC T] = epcmetamodel(constB,constA,k_PL0);

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
DYNAMICS.K_PL0 = k_PL0;
DYNAMICS.PC = PC;
DYNAMICS.LC = LC;
DYNAMICS.T = T;
DYNAMICS.REALTIME = cps;
DYNAMICS.REALDATA = realPC;
DYNAMICS.PCS = pcs;
DYNAMICS.LCS = lcs;
DYNAMICS.ERROR = R;