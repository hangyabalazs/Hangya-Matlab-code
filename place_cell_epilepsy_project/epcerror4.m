function fc = epcerror4(xc)
%EPCERROR4   Error function for EPCFIT4.
%   FC = EPCERROR4(XC) returns least square error in FC. Input argument XC
%   should be 1-by-4 array containing constant 'B', constant 'A', constant
%   'M' and 'k_PL0' parameters for EPCMETAMODEL4.
%
%   FMINSEARCHBND calls EPCERROR4 for minimization.
%
%   See also FMINSEARCHBND, EPCMETAMODEL4 and EPCFIT4.

% Get parameters
constB = xc(1);
constA = xc(2);
constM = xc(3);
k_PL0 = xc(4);

% Call the model
[PC LC T] = epcmetamodel4(constB,constA,constM,k_PL0);

% Count the iterations
global COUNTER
COUNTER = COUNTER + 1;

% Get the series of check points
cps = 60 + [0 30 3*60 3*24*60 3*24*60+30 3*24*60+3*60 7*24*60 7*24*60+30 7*24*60+3*60];
pcs = arrayfun(@(s)PC(find(T<s,1,'last')),cps);
lcs = arrayfun(@(s)LC(find(T<s,1,'last')),cps);

% Real data
realPC = [70.75 53.85 59.80 51.79 39.81 41.67 27.66 13.48 30.12];

% Least square error
R = sum((pcs-realPC).^2);
disp([COUNTER R])

% Model comparison
n = length(realPC);
k = length(xc);
AIC = n * log(R/n) + 2 * k;
AICc = AIC + (2 * k * (k + 1)) / (n - k - 1);
BIC = n * log(R/n) + k * log(n);

% Output
fc = R;
global DYNAMICS
DYNAMICS.CONSTB = constB;
DYNAMICS.CONSTA = constA;
DYNAMICS.CONSTM = constM;
DYNAMICS.K_PL0 = k_PL0;
DYNAMICS.PC = PC;
DYNAMICS.LC = LC;
DYNAMICS.T = T;
DYNAMICS.REALTIME = cps;
DYNAMICS.REALDATA = realPC;
DYNAMICS.PCS = pcs;
DYNAMICS.LCS = lcs;
DYNAMICS.ERROR = R;
DYNAMICS.AIC = AIC;
DYNAMICS.AICC = AICc;
DYNAMICS.BIC = BIC;

global BESTDYNAMICS   % store best fit
if ~isfield(BESTDYNAMICS,'ERROR') || R < BESTDYNAMICS.ERROR
    BESTDYNAMICS.CONSTB = constB;
    BESTDYNAMICS.CONSTA = constA;
    BESTDYNAMICS.CONSTM = constM;
    BESTDYNAMICS.K_PL0 = k_PL0;
    BESTDYNAMICS.PC = PC;
    BESTDYNAMICS.LC = LC;
    BESTDYNAMICS.T = T;
    BESTDYNAMICS.REALTIME = cps;
    BESTDYNAMICS.REALDATA = realPC;
    BESTDYNAMICS.PCS = pcs;
    BESTDYNAMICS.LCS = lcs;
    BESTDYNAMICS.ERROR = R;
    BESTDYNAMICS.AIC = AIC;
    BESTDYNAMICS.AICC = AICc;
    BESTDYNAMICS.BIC = BIC;
end