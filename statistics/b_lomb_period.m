function [pxx,pyy,jmax,prob,z,effm,period_handle] = b_lomb_period(vdisc,bas,lenu)
%LOMB_PERIOD   Lomb-Scargle periodogram from unit after convolution.
%   LOMB_PERIOD(VDISC,BAS,LENU) needs three input arguments: LENU is the length of the original data,
%   VDISC is the discriminate unit and BAS contains the sampling points.
%   
%   VDISC is convolved with Blackmann-Harris function and the convolved data is sampled at the sampling
%   points given in BAS, then LOMBPER is called for the result.
%
%   See also LOMBPER, BLACKMANNHARRIS and PERIODRUN.

% Input arguments check
error(nargchk(3,3,nargin));

% Convolution
z = zeros(1,lenu);
z(vdisc) = 1;
wbh = blackmanharris(1000);
wipunit = conv(z,wbh);

% Lomp periodogram
y = wipunit(bas);
x = bas ./ 10000;
[pxx,pyy,jmax,prob,z,effm,period_handle] = b_lombper(x,y);