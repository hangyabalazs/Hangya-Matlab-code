function imisortframes
%IMISORTFRAMES   Sorts shifted mutual information frames.
%   IMISORTFRAMES loads shifted mutual information frames calculated with
%   different delays (see IMISHIFT), and sorts the frames in the proper
%   order.
%
%   See also IMISHIFT and IMISHIFTCALL.

% Directories
global DATAPATH
eg = '12';
inpdir = [DATAPATH 'Ulbert\EEG_' eg '\MImap\'];
resdir = inpdir;
mm = pwd;

% Load maps and sort frames
load([inpdir 'MIshiftcontrol2_1.mat'])
ld = size(rIMax,3) * 10;
ldb = (size(rIMax,3) - 1) * 10;
chn = size(rIMax,1);
RIM = zeros(chn,chn,ld);
RIML = zeros(chn,chn,ld);
for k = 1:10:ld
    RIM(:,:,k) = rIMax(:,:,(k-1)/10+1);
    RIML(:,:,k) = rIMaxLoc(:,:,(k-1)/10+1);
end
for t = 100:100:900
    load([inpdir 'MIshiftcontrol2_' num2str(t) '.mat'])
    for k = t/100+1:10:ldb
        RIM(:,:,k) = rIMax(:,:,(k-1-t/100)/10+1);
        RIML(:,:,k) = rIMaxLoc(:,:,(k-1-t/100)/10+1);
    end
end

% Save result
rIMax = RIM;
rIMaxLoc = RIML;
cd(resdir)
save MIshiftcontrolfine2 rIMaxLoc rIMax
cd(mm)