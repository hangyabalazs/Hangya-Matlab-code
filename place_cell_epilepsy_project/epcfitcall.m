function [optinitparams optparams optdynamics] = epcfitcall
%EPCFITCALL   Fits kinetic model with different initial parameters.
%   [OPTINITPARAMS OPTPARAMS OPTDYNAMICS] = EPCFITCALL calls EPCFIT with
%   different initial parameters covering the parameter space for linear
%   optimization. It returns the initial parameters which resulted in the
%   best fit (OPTINITPARAMS), parameters of the best fitting model
%   (OPTPARAMS) and a structure containing additional information about the
%   best fitting model (OPTDYNAMICS). It intends to prevent problems from
%   converging to local optima.
%
%   See also EPCFIT and EPCERROR.

% Initial parameters
constB = [0.01 0.03 0.05 0.1 0.5 0.9];
constA = [1 2 5 10 20];
k_PL0 = [0.0002 0.0005 0.001 0.005];

% Fit model
global DYNAMICS
global BESTDYNAMICS
Ri = inf;
for k1 = 1:length(constB)
    for k2 = 1:length(constA)
        for k3 = 1:length(k_PL0)
            disp([k1 k2 k3])
            [cB cA kPL0 Rc] = epcfit(constB(k1),constA(k2),k_PL0(k3));
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