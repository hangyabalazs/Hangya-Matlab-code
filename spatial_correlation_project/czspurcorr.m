function [Rmods CIs] = czspurcorr

% Caller - Load previous data
global DATAPATH     % read excel files
inpdir_xls = [DATAPATH 'Czurko\discriminated2\'];
xlsname = [inpdir_xls 'placedata_new.xls'];
headerrows = 1;
[ntz mtz atz] = xlsread(xlsname,'int');
mtz(1:headerrows,:) = [];
atz(1:headerrows,:) = [];
sf = size(mtz,1);   % number of interneurons
[ntz_pyr mtz_pyr atz_pyr] = xlsread(xlsname,'pyr');
mtz_pyr(1:headerrows,:) = [];
atz_pyr(1:headerrows,:) = [];
Rmods = [];
CIs = [];
for k = 1:sf    % interneuron loop
    animal = atz{k,1};  % rat ID
    
    pyrinx = find(~strcmp(animal,atz_pyr(:,1)));
    sf_pyr = size(pyrinx,1);   % number of pyramidal neurons
    for pk_pyr = 1:sf_pyr    % place cell loop
        k_pyr = pyrinx(pk_pyr);
        animal_pyr = atz_pyr{k_pyr,1};  % rat ID
        if strcmp(animal,animal_pyr)
            error('Technical error 19.')
        end

        [Rmod CI] = czspurcorr_main(atz{k,3},atz_pyr{k_pyr,3},animal,animal_pyr,atz{k,2},atz_pyr{k_pyr,2});
        Rmods(end+1) = Rmod;
        CIs(end+1) = CI;
    end     % end of place cell loop
end     % end of interneuron loop
Rmods = Rmods(~isnan(Rmods));
CIs = CIs(~isnan(CIs));

% -------------------------------------------------------------------------
function [Rmod CI] = czspurcorr_main(int,pyr,exp,exp_pyr,type,type_pyr)

% Load interneuron place map
global DATAPATH
inpdir = [DATAPATH 'Czurko\discriminated2\' type '\' exp '\'];
cd(inpdir)
load('placemaps.mat')
load('placecell_index.mat')
int2 = find(pci==int);
rhst_int = rhst{pci(int2)};

% Load pyramidal cell place map
inpdir = [DATAPATH 'Czurko\discriminated2\' type_pyr '\' exp_pyr '\'];
cd(inpdir)
load('placemaps.mat')
load('placecell_index.mat')
pyr2 = find(pci==pyr);
rhst_pyr = rhst{pci(pyr2)};

% Spatial correlation
Rmod = czpfs_mod(rhst_int,rhst_pyr);
% if Rmod < -0.4
%     keyboard
% end

% Complementarity index
s_int = rhst_int;
ms3_int = s_int .* (zero2nan(double(s_int>b_max_nonnan(s_int)*0.5)));
s_pyr = rhst_pyr;
ms3_pyr = s_pyr .* (zero2nan(double(s_pyr>b_max_nonnan(s_pyr)*0.5)));
c = length(rhst_int(~isnan(rhst_int)));
[C Cb CI] = czcmpl2(ms3_int,ms3_pyr,c);

% -------------------------------------------------------------------------
function R = czpfs_mod(rm1,rm2)
%CZPFS   Place Field Similarity.
%   R = CZPFS_MOD(RM1,RM2) calculates linear correlation between rate maps 
%   RM1 and RM2. It discards areas where neither of the cells fire.
%
%   See also CZPLACEANALYSIS and CZPFS.

% Adjust place filed sizes
[sx1 sy1] = size(rm1);
[sx2 sy2] = size(rm2);
if max(abs(sx1-sx2),abs(sy1-sy2)) >= 5
    R = NaN;
    return
end
sx = min(sx1,sx2);
sy = min(sy1,sy2);
rm1 = rm1(1:sx,1:sy);
rm2 = rm2(1:sx,1:sy);

% Find non-NaNs
r1 = rm1(~isnan(rm1)&~isnan(rm2));
r2 = rm2(~isnan(rm1)&~isnan(rm2));

% Linear correlation
rr1 = r1(r1>b_max_nonnan(r1)*0.2|r2>b_max_nonnan(r2)*0.2);
rr2 = r2(r1>b_max_nonnan(r1)*0.2|r2>b_max_nonnan(r2)*0.2);
corrmtx = corrcoef(rr1,rr2);
R = corrmtx(2);

% -------------------------------------------------------------------------
function [C Cb Cc] = czcmpl2(P,Q,c)
%CZCMPL2   Place field complementarity index.
%   [C CB CC] = CZCMPL(P,Q,L) calculates complementarity between place fields
%   (not rate maps! - see CZPLACEANALYSIS) P (interneuron) and Q (place
%   cell), when L is the size of the arena in pixels. Complementarity index
%   is defined as follows.
%       C = (|Q \ P| / |Q| + |P u Q| / L) / 2
%       CB = (|P u Q|) / (|P| + |Q|) + |P u Q| / L) / 2
%       CC = (|Q \ P| / |Q| + |Q \ P| / (L - |P|)) / 2
%   The first element of the addidion quantifies 'disjunctness', the second
%   stands for coverage.
%
%   See also CZPLACEANALYSIS, CZCMPL and CZPLACEFIGS.

% Adjust place filed sizes
[sx1 sy1] = size(P);
[sx2 sy2] = size(Q);
if max(abs(sx1-sx2),abs(sy1-sy2)) >= 5
    C = NaN;
    Cb = NaN;
    Cc = NaN;
    return
end
sx = min(sx1,sx2);
sy = min(sy1,sy2);
P = P(1:sx,1:sy);
Q = Q(1:sx,1:sy);

% Input argument check
error(nargchk(3,3,nargin))
if ~isequal(size(P),size(Q))
    error('Input distributions must be of the same size.')
end

% Complementarity index
C1 = length(find(isnan(P)&~isnan(Q))) / length(find(~isnan(Q)));     % disjunct
C1b = length(find(~isnan(P)|~isnan(Q))) / (length(find(~isnan(Q))) + length(find(~isnan(P))));
C2 = length(find(~isnan(P)|~isnan(Q))) / c;         % coverage
C2b = length(find(isnan(P)&~isnan(Q))) / (c - length(find(~isnan(P))));
C = (C1 + C2) / 2;
Cb = (C1b + C2) / 2;
Cc = (C1 + C2b) / 2;