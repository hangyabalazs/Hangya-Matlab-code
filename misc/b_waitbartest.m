function b_waitbartest
%WAITBARTEST    Function for experiments on waitbar function.
%
%   See also WAITBAR.

%Input arguments check
error(nargchk(0,0,nargin));

%Waitbar test
h = waitbar(0,'Please wait...','Position',[360 250 275 50]);
for i=1:100,
    pause(.06);
waitbar(i/100)
end
close(h)