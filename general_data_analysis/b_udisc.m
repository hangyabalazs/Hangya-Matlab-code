function [vdisc kuszob] = b_udisc(unit)
%UDISC   Discriminator for raw unit data.
%   [VDISC THR] = UDISC(UNIT) converts UNIT into a series of time values
%   (VDISC) by thresolding (THR).
%
%   See also IN, DISC and DISC2.

% Input arguments check
error(nargchk(1,1,nargin));
unit = unit(:)';    % convert to row vector

% Discriminating
s = figure;
plot(unit,'m');
title('Give the threshold! /Command window/')
kuszob = input('Give the threshold! ');

disc = find(unit>=kuszob);
discl = length(disc);
disc = [disc; unit(disc(1:discl))];
disc(1,discl+1) = length(unit) + 2;
dif = diff(disc(1,:));
difn1 = find(dif>1);
difn1(2:length(difn1)+1) = difn1;
difn1(1) = 0;
vdisc = zeros(1,length(difn1)-1);
for j=1:length(difn1) - 1
    [maxe,maxh] = max(disc(2,difn1(j)+1:difn1(j+1)));
    vdisc(j) = disc(1,difn1(j)+maxh);
end