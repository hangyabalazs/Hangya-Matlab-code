function smooth = gausskernel(sz, arg)
%GAUSSKERNEL  Gaussian kernel for 2D convolution.
%   W = GAUSSKERNEL(SIZE,ARG) sets the size of the convolution kernel
%       (default is [3 3]). If SIZE is a scalar, the size is interpreted
%       as [SIZE SIZE]. ARG sets the standard deviation of the gaussian 
%       kernel (default is .65).
%   
%   See SMOOTH2, SMOOTH3 and SMOOTH2_NONNAN.

% Argument check
if nargin==0
  sz = 3;
  arg = .65;
elseif nargin==1
  arg = .65;
elseif nargin>2
  error(id('WrongNumberOfInputs'),'Wrong number of input arguments.'); 
end

if length(sz)==1
  sz = [sz sz];
elseif numel(sz)~=2
  error(id('InvalidSizeInput'),'SIZE must be a scalar or a 2 element vector.')
end

sz = sz(:)';

padSize = (sz-1)/2;

if ~isequal(padSize, floor(padSize)) || any(padSize<0)
  error(id('InvalidSizeValues'),'All elements of SIZE must be odd integers greater than or equal to 1.');
end

% Kernel
smooth = gaussian2(sz,arg);

% ------------------------------------------------------------------------
function h = gaussian2(P1, P2)
%2D Gaussian lowpass filter
%
%   H = gaussian2(N,SIGMA) returns a rotationally
%   symmetric 2D Gaussian lowpass filter with standard deviation
%   SIGMA (in pixels). N is a 1-by-2 vector specifying the number
%   of rows, columns, pages in H. (N can also be a scalar, in 
%   which case H is NxN.) If you do not specify the parameters,
%   the default values of [3 3] for N and 0.65 for
%   SIGMA.


if nargin>0,
  if ~(all(size(P1)==[1 1]) || all(size(P1)==[1 2])),
     error(id('InvalidFirstInput'),'The first parameter must be a scalar or a 1-by-3 size vector.');
  end
  if length(P1)==1, siz = [P1 P1]; else siz = P1; end
end

if nargin<1, siz = [3 3]; end
if nargin<2, std = .65; else std = P2; end
[x,y] = meshgrid(-(siz(2)-1)/2:(siz(2)-1)/2, -(siz(1)-1)/2:(siz(1)-1)/2);
h = exp(-(x.*x + y.*y)/(2*std*std));
h = h/sum(h(:));

% -------------------------------------------------------------------------
function str=id(str)
str = ['MATLAB:smooth3:' str];
