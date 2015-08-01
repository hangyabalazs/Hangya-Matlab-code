function [optinitparams optparams optdynamics] = epcfitcall6
%EPCFITCALL6   Fits kinetic model with different initial parameters.
%   [OPTINITPARAMS OPTPARAMS OPTDYNAMICS] = EPCFITCALL6 calls EPCFIT6 with
%   different initial parameters covering the parameter space for linear
%   optimization. It returns the initial parameters which resulted in the
%   best fit (OPTINITPARAMS), parameters of the best fitting model
%   (OPTPARAMS) and a structure containing additional information about the
%   best fitting model (OPTDYNAMICS). It intends to prevent problems from
%   converging to local optima.
%
%   See also EPCFIT6 and EPCERROR6.

% Initial parameters
constB = [0.006 0.06 0.6];
constA = [1 7 21];
constM = [15 50];
k_LP0 = [0.00001 0.0002 0.005];

% Fit model
global DYNAMICS
global BESTDYNAMICS
Ri = inf;
for k1 = 1:length(constB)
    for k2 = 1:length(constA)
        for k3 = 1:length(constM)
            for k4 = 1:length(k_LP0)
                disp([k1 k2 k3 k4])
                [cB cA cM kLP0 Rc] = epcfit6(constB(k1),constA(k2),constM(k3),k_LP0(k3));
                if Rc < Ri
                    Ri = Rc;
                    optinitparams = [constB(k1) constA(k2) constM(k3) k_LP0(k3)];
                    optparams = [cB cA cM kLP0];
                    optdynamics = DYNAMICS;
                    save tmp optinitparams optparams optdynamics DYNAMICS BESTDYNAMICS k1 k2 k3 k4
                end
            end
        end
    end
end