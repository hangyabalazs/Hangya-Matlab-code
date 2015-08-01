function snwfreqvscoupling(allsegments)
%SNWFREQVSCOUPLING   Frequency of coupling types.
%   SNWFREQVSCOUPLING calculates the relative frequency of different types
%   of sniff-whisk phase locking (see SNWPHASE4 and SNWSNIFFPLOT) and plots
%   them against sniffing frequency.
%
%   SNWFREQVSCOUPLING(ALLSEGMENTS) uses the segments structure given in the
%   input argument.
%
%   See also SNWPHASE4 and SNWSNIFFPLOT.

% Input argument check
error(nargchk(0,1,nargin))
if nargin < 1
    load('C:\Balazs\_analysis\SniffWhisk\Fig_frequency\allsegments.mat')
end

% Frequency of coupling types - pool data
L = zeros(5,15);    % dimensions: coupling type, sniff freq.
for k = 1:length(allsegments)
    cfr = max(1,min(floor(allsegments(k).frequency),15));   % cut sniff freq at 1 and 15
    switch allsegments(k).group
        case 1/3
            L(1,cfr) = L(1,cfr) + allsegments(k).length;
        case 0.5
            L(2,cfr) = L(2,cfr) + allsegments(k).length;
        case 1
            L(3,cfr) = L(3,cfr) + allsegments(k).length;
        case 2
            L(4,cfr) = L(4,cfr) + allsegments(k).length;
        case 3
            L(5,cfr) = L(5,cfr) + allsegments(k).length;
    end
end
Lnorm = L / sum(L(3,:));    % normalize to 1:1 coupling
figure;plot(Lnorm')

% Frequency of coupling types - average for animals first
L = zeros(5,15,4);    % dimensions: coupling type, sniff freq., animal
for k = 1:length(allsegments)
    cfr = max(1,min(floor(allsegments(k).frequency),15));   % cut sniff freq at 1 and 15
    switch allsegments(k).rat
        case 'R1'
            ratinx = 1;
        case 'R3'
            ratinx = 2;
        case 'R7'
            ratinx = 3;
        case 'P5'
            ratinx = 4;
    end
            
    switch allsegments(k).group
        case 1/3
            L(1,cfr,ratinx) = L(1,cfr,ratinx) + allsegments(k).length;
        case 0.5
            L(2,cfr,ratinx) = L(2,cfr,ratinx) + allsegments(k).length;
        case 1
            L(3,cfr,ratinx) = L(3,cfr,ratinx) + allsegments(k).length;
        case 2
            L(4,cfr,ratinx) = L(4,cfr,ratinx) + allsegments(k).length;
        case 3
            L(5,cfr,ratinx) = L(5,cfr,ratinx) + allsegments(k).length;
    end
end
pLnorm = zeros(size(L));
for k = 1:4    % normalize to 1:1 coupling, animal by animal
    pLnorm(:,:,k) = L(:,:,k) / sum(L(3,:,k));
end
Lnorm = squeeze(mean(pLnorm,3));
figure
plot(Lnorm')