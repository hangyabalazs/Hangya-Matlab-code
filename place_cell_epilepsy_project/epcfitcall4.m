function [optinitparams optparams optdynamics] = epcfitcall4
%EPCFITCALL4   Fits kinetic model with different initial parameters.
%   [OPTINITPARAMS OPTPARAMS OPTDYNAMICS] = EPCFITCALL4 calls EPCFIT4 with
%   different initial parameters covering the parameter space for linear
%   optimization. It returns the initial parameters which resulted in the
%   best fit (OPTINITPARAMS), parameters of the best fitting model
%   (OPTPARAMS) and a structure containing additional information about the
%   best fitting model (OPTDYNAMICS). It intends to prevent problems from
%   converging to local optima.
%
%   See also EPCFIT4 and EPCERROR4.

% Initial parameters
constB = [0.002 0.05 0.7];
constA = [0.6 1 20];
constM = [0.1 0.7];
k_PL0 = [0.0002 0.005 0.07];

% Fit model
global DYNAMICS
global BESTDYNAMICS
Ri = inf;
for k1 = 1:length(constB)
    for k2 = 1:length(constA)
        for k3 = 1:length(constM)
            for k4 = 1:length(k_PL0)
                disp([k1 k2 k3 k4])
                [cB cA cM kPL0 Rc] = epcfit4(constB(k1),constA(k2),constM(k3),k_PL0(k4));
                if Rc < Ri
                    Ri = Rc;
                    optinitparams = [constB(k1) constA(k2) constM(k3) k_PL0(k4)];
                    optparams = [cB cA cM kPL0];
                    optdynamics = DYNAMICS;
                    save tmp optinitparams optparams optdynamics DYNAMICS BESTDYNAMICS k1 k2 k3 k4
                end
            end
        end
    end
end