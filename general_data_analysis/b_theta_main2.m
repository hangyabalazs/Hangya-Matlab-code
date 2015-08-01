function b_theta_main2
%THETA_MAIN2    Ultimate theta selector.
%   THETA_MAIN2 calls the thetaselector functions. It performs the entire theta
%   segment selection.
%
%   THETA_MAIN2 works on two directories: the raw data directory and the results'
%   directory. You are able to modify them through editing the program code.
%
%   See also THETASELECTORRUN2, THETA and NONTHETA.

% Input argument check
error(nargchk(0,0,nargin));

% Define directories
global DATAPATH     % b_startup is required!
global DATADIR
inpdir = [DATADIR 'raphe_matfiles\raphe_juxta_files_discriminated\temp\'];
resdir = [DATAPATH 'raphe\raphe_juxta\Wavelet_temp\'];
cd(resdir)

% Thetaselectorrun
disp('Running THETASELECTORRUN2...');
b_thetaselectorrun2(inpdir,resdir);

% Load scalevector and newstep
load scalevector.mat
load newstep.mat

% Theta
disp('Running THETA...');
b_theta(f,newstep,resdir)

% Nontheta
disp('Running NONTHETA...')
b_nontheta(f,newstep,resdir,inpdir)