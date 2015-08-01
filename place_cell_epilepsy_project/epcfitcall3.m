function [optinitparams optparams optdynamics] = epcfitcall3
%EPCFITCALL3   Fits kinetic model with different initial parameters.
%   [OPTINITPARAMS OPTPARAMS OPTDYNAMICS] = EPCFITCALL3 calls EPCFIT3 with
%   different initial parameters covering the parameter space for linear
%   optimization. It returns the initial parameters which resulted in the
%   best fit (OPTINITPARAMS), parameters of the best fitting model
%   (OPTPARAMS) and a structure containing additional information about the
%   best fitting model (OPTDYNAMICS). It intends to prevent problems from
%   converging to local optima.
%
%   See also EPCFIT3 and EPCERROR3.

% Initial parameters
constB = [0.05 5 20];
constA = [0.5 1 10 50];
k_PL0 = [0.0002 0.005 0.07];

% Fit model
global DYNAMICS
global BESTDYNAMICS
Ri = inf;
for k1 = 1:length(constB)
    for k2 = 1:length(constA)
        for k3 = 1:length(k_PL0)
            disp([k1 k2 k3])
            [cB cA kPL0 Rc] = epcfit3(constB(k1),constA(k2),k_PL0(k3));
            if Rc < Ri
                Ri = Rc;
                optinitparams = [constB(k1) constA(k2) k_PL0(k3)];
                optparams = [cB cA kPL0];
                optdynamics = DYNAMICS;
                save tmp optinitparams optparams optdynamics DYNAMICS BESTDYNAMICS k1 k2 k3
            end
        end
    end
end