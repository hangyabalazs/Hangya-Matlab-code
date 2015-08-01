function b_theta_main
%THETA_MAIN    Ultimate theta selector.
%   THETA_MAIN calls the thetaselector functions. It performs the entire theta
%   segment selection.
%
%   THETA_MAIN works on two directories: the raw data directory and the results'
%   directory. You are able to modify them through editing the program code.
%
%   See also THETASELECTORRUN, THETASELECTORRUN_LONG, LONGCOPY, THETA, THETA_LONG,
%   NONTHETA and NONTHETA_LONG.

% Input argument check
error(nargchk(0,0,nargin));

% Define directories
global DATAPATH     % b_startup is required!
global DATADIR
inpdir = [DATAPATH 'DATA\analysenow5\'];
resdir = [DATAPATH 'Wavelet2\'];
cd(resdir)

% Thetaselectorrun
disp('Running THETASELECTORRUN...');
long = b_thetaselectorrun(inpdir,resdir);

% Longcopy
disp('Running LONGCOPY...');
b_longcopy(long,inpdir)

% Thetaselectorrun_long
disp('Running THETASELECTORRUN_LONG...');
[f,newstep] = b_thetaselectorrun_long(inpdir,resdir);

% Theta
disp('Running THETA...');
b_theta(f,newstep,resdir)

% Theta_long
disp('Running THETA_LONG...')
b_theta_long(f,newstep,resdir)

% Nontheta
disp('Running NONTHETA...')
b_nontheta(f,newstep,resdir,inpdir)

% Nontheta_long
disp('Running NONTHETA_LONG...')
b_nontheta_long(f,newstep,resdir,inpdir)