function sl = imisiglev(cI)
%IMISIGLEV   Significance levels for MI calculations.
%   SL = IMISIGLEV(CI) calculates critical values for different
%   significance levels for mutual information distributions using a
%   control set given in CI.
%
%   To get significance levels for IMISHIFT result files, CI should be
%   calculated using IMIRAW and IMIRAWCALL.
%
%   See also IMIRAWCALL.

% Input argumnet check
error(nargchk(0,1,nargin))
ls = 0;     % indicator for load and save
if nargin < 1
    global DATAPATH
    eg = '53';
    dr = [DATAPATH 'Ulbert\OITI_37_EEG_' eg '\CCGmap\'];
    cd(dr)
    load([dr 'MIrawmaps_EEG' eg '_filt01_40_overlap10.mat'])
    ls = 1;
end

% Critical values
for x = 1:20;
    for y = 1:20
        cI = cI(:);
        cI = cI(cI>0);
        sci = sort(cI(:));
        lci = length(cI(:));
        levs = [0.95 0.99 0.999 0.9999];  % significance levels
        lle = length(levs);
        sl = zeros(lle,2);
        for k = 1:lle
            sl(k,1) = levs(k);
            inx = ceil(lci*levs(k));
            sl(k,2) = sci(inx);     % critical values
        end
    end
end

% Save
if ls
    fn = ['siglev_EEG' eg '.mat'];
    save(fn,'sl')
end