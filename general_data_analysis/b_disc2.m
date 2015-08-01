function b_disc2(kuszob2,seglen)
%DISC2   Discriminator for raw unit data.
%   DISC2 converts unit into a series of time values.
%   The output of DISC is a global variable with the same name. This variable is a cell array containing
%   1. id: the matlab id of the input raw data (id{1}:name, id{2}:fisrt point of interval, id{3}:last 
%   point of interval),
%   2. output: 1st row contains the time series, 2nd row contains the value of raw data at the maximum, 
%   of spike; the last element of the matrix is the end of discriminated part,
%   3. vdisc: contains only the time series,
%   4. kuszob: the limit amplitude value for discrimination,
%   5. instfrek: instant frequency,
%   6. isi: interspike intervalls (ordered in time).
%
%   DISC2(KUSZOB) does the discrimination using KUSZOB as threshold.
%
%   DISC2, unlike DISC can deal with 'changing threshold' created by THRES2. Use DISC2(KUSZOB,SEGLEN)
%   syntax in this case where SEGLEN is the segment length by thresholding - the second output argument
%   of THRES2.
%
%   Note, that DISC works only if the data is inported!
%
%   See also IN, THRES and THRES2.

% Input arguments check
error(nargchk(0,2,nargin));

% Reading the input data
global IN
if isempty(IN)
    error('Data has not been inported.')
end
data = IN{1};
eeg = IN{2};
fname = IN{3};
pathname = IN{4};
datinx1 = IN{5};
datinx2 = IN{6};
time = IN{7};
unit = IN{8};
dt = IN{9};
meret = IN{10};
mintafr = IN{11};
xlimit = IN{12};

% Check if the datafile is already discriminated
global DISC
if isempty(DISC) == 0,
    id = DISC{1};
    output = DISC{2};
    vdisc = DISC{3};
    kuszob = DISC{4};
    instfrek = DISC{5};
    isi = DISC{6};
    if nargin == 0,
        if isequal(id{1},fname) & isequal(id{2},datinx1) & isequal(id{3},datinx2),
            return;
        end;
    end;
    if nargin == 1,
        if isequal(id{1},fname) & isequal(id{2},datinx1) & isequal(id{3},datinx2) & isequal(kuszob,kuszob2),
            return;
        end;
    end;
    clear global DISC
end;

% Threshold
switch nargin
case 0
    s = figure;
    plot(unit,'m');
    title('Give the threshold! /Command window/')
    kuszob = input('Give the threshold! ');
case 1
    kuszob = kuszob2;
case 2
    ind2 = 0;
    lenu = length(unit);
    next = 1;
    kuszob = [];
    while ind2 < lenu
        ind1 = ind2 + 1;
        ind2 = ind2 + seglen;
        kuszob(ind1:min(ind2,lenu)) = kuszob2(next);
        next = next + 1;
    end
end

% Discriminating
disc = find(unit>=kuszob); 
discl = length(disc);
disc = [disc; unit(disc(1:discl))];
disc(1,discl+1) = length(unit) + 2;
dif = diff(disc(1,:));
difn1 = find(dif>1);
difn1(2:length(difn1)+1) = difn1;
difn1(1) = 0;
vdisc = zeros(1,length(difn1)-1);
for j = 1:length(difn1) - 1
    [maxe,maxh] = max(disc(2,difn1(j)+1:difn1(j+1)));
    vdisc(j) = disc(1,difn1(j)+maxh);
end

% Computing the instant frequency
if ~vdisc == 0,
    isi = diff(vdisc) * dt;
    instfrek = zeros(1,length(unit));
    for i=1:length(vdisc)-1
        instfrek(vdisc(i):vdisc(i+1)) = 1 / isi(i);
    end;
    instfrek(1:vdisc(1)-1) = 1 / (vdisc(1) * dt);
    instfrek(vdisc(end):length(unit)) = 1 / ((length(unit) - vdisc(end)) * dt);
%    bins = 0.01:0.01:1;
%    [isih,xout] = hist(isi,bins);
        
    output = [vdisc;unit(vdisc(:))];
    if output(1,end) ~= length(unit), 
        output(1,end+1) = length(unit); 
    end;
else output = length(unit);
end
if exist('s'),
    close(s);
end;

% Creating id
id = cell(1,3);
id{1} = fname;
id{2} = datinx1;
id{3} = datinx2;

% Creating the global variable DISC
global DISC
DISC = cell(1,6);
DISC{1} = id;
DISC{2} = output;
DISC{3} = vdisc;
DISC{4} = kuszob;
DISC{5} = instfrek;
DISC{6} = isi;