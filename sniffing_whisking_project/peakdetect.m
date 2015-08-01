function [ind, vpeaks] = peakdetect(V,threshold)
%PEAKDETECT   Find peaks of a signal.
%  [ind, peaks] = peakdetect(Signal,threshold);
%
% AK 2000/1

V=V(:);

dv=diff(V);                 %first derivative
ddv=diff(dv);               %second derivative
dv1=[dv(2:length(dv))' 0]'; %shift if over

w=dv.*dv1;                  %look at neighbour's behavior
pw=find(w<0);               %if sign change occured we're at home	

ps=pw(find(ddv(pw)<0));             %choose only the max based on 2nd der test
ind=ps(find(V(ps+1)>threshold)) +1;	%threshold & correct spike timing 
vpeaks = V(ind);
