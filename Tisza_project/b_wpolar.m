function b_wpolar(ang)
%WPOLAR   Plots polar plot, rose histogram and compass from phase values.
%   WPOLAR(ANG) plots polar plot, angle histogram and compass from phase
%   values in ANG.
%
%   See also WANGLE2.

% Polar plot
n = length(ang);
figure
H = polar(ang,ones(1,n),' ro');
set(H,'MarkerSize',5)

% Rose diagram
figure
rose(ang)

% Compass
figure
compass(exp(1).^(i.*ang))