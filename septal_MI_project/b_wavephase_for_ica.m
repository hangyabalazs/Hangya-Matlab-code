function [th_int,no_int] = wavephase_for_ica(wave,f,newstep)
%WAVEPHASE_FOR_ICA Version of wavephase used by ICA_WAVE.
%   [TH_INT, NO_INT] = WAVEPHASE_FOR_ICA(WAVE,F,NEWSTEP) needs three input arguments:
%   WAVELET output (WAVE), scalevector (F) and downsampling constant (NEWSTEP) for 10
%   kHz sampled data. It returns theta (TH_INT) and non-theta (NO_INT) intervals.
%
%   WAVEPHASE_FOR_ICA uses the following theta definition: more than 30 % of the eeg
%   wavelet power fall in the 3 - 6 Hz band.
%
%   See also WAVEPHASE and ICA.

% Input arguments check
error(nargchk(3,3,nargin));

% Input variables
global IN
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

% Wavelet power calculation
n = fix(length(eeg)/newstep) - 1;
power = (abs(wave)) .^ 2;
fnd = find(f>6);
pwind1 = fnd(end); %if pwind1 ~= 76, error('pwind1 error'); end
fnd = find(f<3);
pwind2 = fnd(1); %if pwind2 ~= 102, error('pwind2 error'); end

% Selection of segments
for segment_type = 1:2
    thp = zeros(1,n);
    allp = zeros(1,n);
    selector = zeros(1,n);
    for i = 1:n
        thp(1,i) = sum(power(pwind1:pwind2,i));
        allp(1,i) = sum(power(:,i));
    end
    thprop = thp ./ allp;
    switch segment_type
    case 1  % Selection of segments w. >30% ~3-6 Hz frequency content
        for i = 1:n
            if thprop(1,i) > 0.3
                selector(1,i) = 1;
            end
        end
    case 2  % Selection of segments w. <30% ~3-6 Hz frequency content
        for i = 1:n
            if thprop(1,i) < 0.3
                selector(1,i) = 1;
            end
        end
    end
    ds = diff(selector);
    ds1 = zeros(1,n-1);
    ds1 = abs(ds);
    if selector(1) == 1
        ds1(1) = 1;
    end
    if selector(end) == 1
        ds1(end) = 1;
    end
    ds2 = find(ds1(1,:)==1);
    if rem(length(ds2),2) == 1
        error('Length ds2 is odd!');
    end
    
    e1 = diff(ds2);   %leaving the intervals shorter than 5 sec...
    e2 = e1(1:2:end);
    e3 = find(e2>50000/newstep);
    e4 = zeros(1,2*length(e3));
    for i = 1:length(e3)
        e4(2*i) = 2*e3(i);
        e4(2*i-1) = e4(2*i) - 1;
    end
    e5 = ds2(e4);
    ds3 = e5;
    
    ds4 = ds3 .* newstep + datinx1;    %localization of first and last points of theta intervals
    ds5 = ds4(1,1:2:end); %first points of theta intervals
    ds6 = ds4(1,2:2:end); %last points of theta intervals
    switch segment_type
    case 1
        th_int = ds4;
    case 2
        no_int = ds4;
    end
end