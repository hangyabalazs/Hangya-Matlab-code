function mat2dat
%MAT2DAT   Converts mat files into dat files.
%   MAT2DAT opens a mat file containing multichannel EEG data and saves
%   data content into a dat file, concatenating the data from all channels.
%   An nfo file is also saved containing the number of channels and the
%   length of the data.

% Import
fn = 'F:\raw_data\human_SO\oiti37_lukacs\grid\mat_EEG_12b\EEG_12_1855_1915_rs_filt01_40.mat';
[pthname fname ext] = fileparts(fn);
datname = [pthname '\' fname '.dat'];
infdat = [pthname '\' fname '.nfo'];
load(fn)

% Write dat
chno = size(data,2);    % number of channels
datalen = size(data,1);    % length of data
fid = fopen(datname,'w');
for k = 1:chno
    fwrite(fid,data(:,k),'double');
end
fclose(fid);

% Write nfo
fid = fopen(infdat,'w');
fwrite(fid,chno,'double');
fwrite(fid,datalen,'double');
fclose(fid);