function msptheta_main2
%MSPTHETA_MAIN2    Ultimate theta selector.
%   MSPTHETA_MAIN2 is modified from THETA_MAIN2 (see the help of the latter
%   for details).
%
%   See also MSPTHETASELECTORRUN2, THETA and NONTHETA.

% Input argument check
error(nargchk(0,0,nargin));

% Define directories
global DATAPATH     % b_startup is required!
global DATADIR
inpdir = [DATADIR 'MSHCsp\viktor5\EEG\'];
resdir = [DATAPATH 'MSHCsp\Wavelet_temp\'];
cd(resdir)

% Samplng rate
sr = 20000;
disp(['Sampling rate: ' num2str(sr)])

% Thetaselectorrun
disp('Running MSPTHETASELECTORRUN2...');
mspthetaselectorrun2(inpdir,resdir,sr);

% Load scalevector and newstep
load scalevector.mat
load newstep.mat

% Theta
disp('Running MSPTHETA...');
msptheta(f,newstep,resdir)