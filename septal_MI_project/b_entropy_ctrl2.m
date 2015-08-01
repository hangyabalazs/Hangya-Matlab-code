function b_entropy_ctrl2
%ENTROPY_CTRL2   Controls for entropy project.
%   ENTROPY_CTRL2 performs 6 type of entropy calculation calling ENTROPY function:
%       1. random data vs. magic square
%       2. random data vs. constant function
%       3. random data vs. linear function
%       4. random data vs. polinome function
%       5. random data vs. exponential function
%       6. random data vs. itself + 1.
%
%   ENTROPY_CTRL2 calulates entropy on larger sets than ENTROPY_CTRL.
%
%   See also ENTROPY and ENTROPY_CTRL.

% Input argument check
error(nargchk(0,0,nargin))

% Define directories
global DATAPATH
resdir = [DATAPATH 'Entropy\control\control_10000\'];      % results' directory
mm = pwd;
cd(resdir)
create_subdir

% 'RANDOM' vs. 'MAGIC'
RHx = []; RHy = []; RHxy = [];
RHxcy = []; RHycx = [];
RIxy = []; RUxy = []; RUyx = [];
for i = 1:100
    y1 = rand(1,10000);
    y2 = magic(100);
    y2 = y2(:)';
    [hx,hy,jh,Hx,Hy,Hxy,Hxcy,Hycx,Ixy,Uxy,Uyx] = b_entropy(y1,y2);
    
    RHx = [RHx Hx];
    RHy = [RHy Hy];
    RHxy = [RHxy Hxy];
    RHxcy = [RHxcy Hxcy];
    RHycx = [RHycx Hycx];
    RIxy = [RIxy Ixy];
    RUxy = [RUxy Uxy];
    RUyx = [RUyx Uyx];
end
fln = 'RANDOM_MAGIC';
cd random_magic
save(fln,'RHx','RHy','RHxy','RHxcy','RHycx','RIxy','RUxy','RUyx')
cd ..

% 'RANDOM' vs. 'CONSTANT'
RHx = []; RHy = []; RHxy = [];
RHxcy = []; RHycx = [];
RIxy = []; RUxy = []; RUyx = [];
for i = 1:100
    y1 = rand(1,10000);
    y2 = ones(1,10000);
    y2 = y2(:)';
    [hx,hy,jh,Hx,Hy,Hxy,Hxcy,Hycx,Ixy,Uxy,Uyx] = entropy(y1,y2);
    
    RHx = [RHx Hx];
    RHy = [RHy Hy];
    RHxy = [RHxy Hxy];
    RHxcy = [RHxcy Hxcy];
    RHycx = [RHycx Hycx];
    RIxy = [RIxy Ixy];
    RUxy = [RUxy Uxy];
    RUyx = [RUyx Uyx];
end
fln = 'RANDOM_CONSTANT';
cd random_constant
save(fln,'RHx','RHy','RHxy','RHxcy','RHycx','RIxy','RUxy','RUyx')
cd ..

% 'RANDOM' vs. 'LINEAR'
RHx = []; RHy = []; RHxy = [];
RHxcy = []; RHycx = [];
RIxy = []; RUxy = []; RUyx = [];
for i = 1:100
    y1 = rand(1,10000);
    y2 = [1:10000];
    y2 = y2(:)';
    [hx,hy,jh,Hx,Hy,Hxy,Hxcy,Hycx,Ixy,Uxy,Uyx] = b_entropy(y1,y2);
    
    RHx = [RHx Hx];
    RHy = [RHy Hy];
    RHxy = [RHxy Hxy];
    RHxcy = [RHxcy Hxcy];
    RHycx = [RHycx Hycx];
    RIxy = [RIxy Ixy];
    RUxy = [RUxy Uxy];
    RUyx = [RUyx Uyx];
end
fln = 'RANDOM_LINEAR';
cd random_linear
save(fln,'RHx','RHy','RHxy','RHxcy','RHycx','RIxy','RUxy','RUyx')
cd ..

% 'RANDOM' vs. 'POLINOME'
RHx = []; RHy = []; RHxy = [];
RHxcy = []; RHycx = [];
RIxy = []; RUxy = []; RUyx = [];
for i = 1:100
    y1 = rand(1,10000);
    x =  [1:10000];
    y2 = x .^ 4 + 8 * x .^ 3 - 2 * x + 7;
    y2 = y2(:)';
    [hx,hy,jh,Hx,Hy,Hxy,Hxcy,Hycx,Ixy,Uxy,Uyx] = b_entropy(y1,y2);
    
    RHx = [RHx Hx];
    RHy = [RHy Hy];
    RHxy = [RHxy Hxy];
    RHxcy = [RHxcy Hxcy];
    RHycx = [RHycx Hycx];
    RIxy = [RIxy Ixy];
    RUxy = [RUxy Uxy];
    RUyx = [RUyx Uyx];
end
fln = 'RANDOM_POLINOME';
cd random_polinome
save(fln,'RHx','RHy','RHxy','RHxcy','RHycx','RIxy','RUxy','RUyx')
cd ..

% 'RANDOM' vs. 'EXPONENTIAL'
RHx = []; RHy = []; RHxy = [];
RHxcy = []; RHycx = [];
RIxy = []; RUxy = []; RUyx = [];
for i = 1:100
    y1 = rand(1,10000);
    x =  [1:10000];
    y2 = exp(x/10000);
    y2 = y2(:)';
    [hx,hy,jh,Hx,Hy,Hxy,Hxcy,Hycx,Ixy,Uxy,Uyx] = b_entropy(y1,y2);
    
    RHx = [RHx Hx];
    RHy = [RHy Hy];
    RHxy = [RHxy Hxy];
    RHxcy = [RHxcy Hxcy];
    RHycx = [RHycx Hycx];
    RIxy = [RIxy Ixy];
    RUxy = [RUxy Uxy];
    RUyx = [RUyx Uyx];
end
fln = 'RANDOM_EXPONENTIAL';
cd random_exponential
save(fln,'RHx','RHy','RHxy','RHxcy','RHycx','RIxy','RUxy','RUyx')
cd ..

% 'RANDOM' vs. 'RANDOM + 1'
RHx = []; RHy = []; RHxy = [];
RHxcy = []; RHycx = [];
RIxy = []; RUxy = []; RUyx = [];
for i = 1:100
    y1 = rand(1,1000);
    y2 = y1 + 1;
    [hx,hy,jh,Hx,Hy,Hxy,Hxcy,Hycx,Ixy,Uxy,Uyx] = b_entropy(y1,y2);
    
    RHx = [RHx Hx];
    RHy = [RHy Hy];
    RHxy = [RHxy Hxy];
    RHxcy = [RHxcy Hxcy];
    RHycx = [RHycx Hycx];
    RIxy = [RIxy Ixy];
    RUxy = [RUxy Uxy];
    RUyx = [RUyx Uyx];
end
fln = 'RANDOM_ITSELF';
cd random_itself
save(fln,'RHx','RHy','RHxy','RHxcy','RHycx','RIxy','RUxy','RUyx')
cd ..

cd(mm)



% ----------------------------------------------------------------------------------------
function [hx,hy,jh,Hx,Hy,Hxy,Hxcy,Hycx,Ixy,Uxy,Uyx] = entropy(y1,y2)
%ENTROPY   Entropy for 'RANDOM' vs. 'CONSTANT'.
%   [hx,hy,jh,Hx,Hy,Hxy,Hxcy,Hycx,Ixy,Uxy,Uyx] = ENTROPY(y1,y2) needs 2
%   inputs:
%       y1: first data series
%       y2: second data series
%
%   It produces 11 outputs:
%       hx: normalized distribution of y1
%       hy: normalized distribution of y2
%       jh: common distribution of y1 and y2
%       Hx: y1 wavelet magnitude entropy
%       Hy: y2 wavelet magnitude entropy
%       Hxy: combined entropy
%       Hxcy: conditional entropy (H(y1|y2))
%       Hycx: conditional entropy (H(y2|y1))
%       Ixy: mutual information
%       Uxy: uncertainity coefficient (y2->y1)
%       Uyx: uncertainity coefficient (y1->y2)

% Input argument check
error(nargchk(2,2,nargin))

% Histogram estimation
n1 = size(y1,1) * size(y1,2);    % size of data
n2 = size(y2,1) * size(y2,2);

miny1 = min(min(y1));     % minimum of data1
maxy1 = max(max(y1));     % maximum of data1
h1 = fix(exp(0.626+0.4*log(n1-1)));   % number of bins
binwidth1 = (maxy1 - miny1) ./ h1;
xx1 = miny1 + binwidth1 * (0:h1);   % bin edges
xx1(length(xx1)) = maxy1;
xx1(1) = -inf;
x1 = xx1(1:length(xx1)-1) + binwidth1 / 2;     % bin halves


miny2 = min(min(y2));     % minimum of data2
maxy2 = max(max(y2));     % maximum of data2
h2 = fix(exp(0.626+0.4*log(n2-1)));   % number of bins
binwidth2 = (maxy2 - miny2) ./ h2;
xx2 = miny2 + binwidth2 * (0:h2);   % bin edges
xx2(length(xx2)) = maxy2;
xx2(1) = -inf;
x2 = xx2(1:length(xx2)-1) + binwidth2 / 2;     % bin halves

nbin1 = length(xx1);
nbin2 = length(xx2);
jh = zeros(nbin1-1,nbin2-1);

ty1 = y1 - miny1;
ty2 = y2 - miny2;
p = fix(ty1/binwidth1-10000*eps) + 1;
q = ones(size(y2));
for i = 1:n1
    jh(p(i),q(i)) = jh(p(i),q(i)) + 1;
end

% Calculation of entropies & uncertainity coefficients
[m,n] = size(jh); 
N = sum(sum(jh));
hxx = sum(jh');      % marginal distribution: y1
hyy = sum(jh);     % marginal distribution: y2
N1 = sum(hxx);
N2 = sum(hyy);
hx = hxx / N1;      % normalized distribution: y1
hy = hyy / N2;      % normalized distribution: y2

Hx = 0;
for i = 1:m
    if hx(i) ~= 0
        a = hxx(i) / N;     % normalization
        Hx = Hx - (a * log(a));   % Y1 ENTROPY
    end
end

Hy = 0;
for k = 1:n
    if hy(k) ~= 0
        a = hyy(k) / N;     % normalization
        Hy = Hy - (a * log(a));   % Y2 ENTROPY
    end
end

Hxy = 0;
for i = 1:m
    for k = 1:n 
        if jh(i,k) ~= 0
            a = jh(i,k) / N;     % normalization
            Hxy = Hxy - a * log(a);     % COMMON ENTROPY
        end
    end
end

Hycx = Hxy - Hx;    % conditional entropy
Hxcy = Hxy - Hy;

Uyx = (Hy - Hycx) / (Hy + eps);     % uncertainity coefficient
Uxy = (Hx - Hxcy) / (Hx + eps);
Ux2y = 2 * ((Hy + Hx - Hxy) / (Hx + Hy + eps));

Ixy = Hx + Hy - Hxy;    % mutual information



% ----------------------------------------------------------------------------------------
function create_subdir

if ~b_isdir2('random_magic')
    mkdir random_magic
end
if ~b_isdir2('random_constant')
    mkdir random_constant
end
if ~b_isdir2('random_linear')
    mkdir random_linear
end
if ~b_isdir2('random_polinome')
    mkdir random_polinome
end
if ~b_isdir2('random_exponential')
    mkdir random_exponential
end
if ~b_isdir2('random_itself')
    mkdir random_itself
end