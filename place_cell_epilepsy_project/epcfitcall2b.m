function [optinitparams optparams optdynamics] = epcfitcall2b
%EPCFITCALL2B   Fits kinetic model with different initial parameters.
%   [OPTINITPARAMS OPTPARAMS OPTDYNAMICS] = EPCFITCALL2B calls EPCFIT2B with
%   different initial parameters covering the parameter space for linear
%   optimization. It returns the initial parameters which resulted in the
%   best fit (OPTINITPARAMS), parameters of the best fitting model
%   (OPTPARAMS) and a structure containing additional information about the
%   best fitting model (OPTDYNAMICS). It intends to prevent problems from
%   converging to local optima.
%
%   See also EPCFIT2B and EPCERROR2B.

% Initial parameters
constM = [15 40 70];
k_PL0 = [0.00002 0.0002 0.0005];

% Fit model
global DYNAMICS
global BESTDYNAMICS
Ri = inf;
for k1 = 1:length(constM)
    for k2 = 1:length(k_PL0)
        disp([k1 k2])
        [cM kPL0 Rc] = epcfit2b(constM(k1),k_PL0(k2));
        if Rc < Ri
            Ri = Rc;
            optinitparams = [constM(k1) k_PL0(k2)];
            optparams = [cM kPL0];
            optdynamics = DYNAMICS;
            save tmp2b optinitparams optparams optdynamics DYNAMICS BESTDYNAMICS k1 k2
        end
    end
end