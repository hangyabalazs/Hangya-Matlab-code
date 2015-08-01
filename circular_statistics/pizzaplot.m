function pizzaplot(a1,a2,am,c)
%PIZZAPLOT  Slice on a polar plot.
%   PIZZAPLOT(A1,A2,AM,C) draws a pizza slice shape patch of color C between
%   angles A1 and A2 via AM (in rad). AM determines which of the two
%   complementary slices to color.
%
%   See also POLAR.

%   Balazs Hangya
%   Cold Spring Harbor Laboratory
%   1 Bungtown Road, Cold Spring Harbor
%   New York 11724, USA
%   balazs.cshl@gmail.com
%
%   Institute of Experimental Medicine
%   Szigony street 43, Budapest 1083, Hungary
%   hangyab@koki.hu

% Decide direction
if (am - a1) * (am - a2) > 0   % am is not between a1 and a2
     sa = sort([a1 a2],'ascend');  % assure a1 < a2
     [a1 a2] = deal(sa(1),sa(2));
     if am > a2  % am above the range
         a1 = a1 + 2 * pi;
     elseif am < a1  % am below the range
         a2 = a2 - 2 * pi;
     end
end

% Angle and radial distance values for the plot
incr = sign(a2-a1) * pi/256;
ang = [a1 a1:incr:a2 a2];
rho = [0 ones(1,length(ang)-2) 0];

% Draw the outline to set up the polar plot
P = polar(ang,rho);
set(P,'Color',c)

% Fill in the patch
[pang prho] = pol2cart(ang,rho);
patch(pang,prho,c)