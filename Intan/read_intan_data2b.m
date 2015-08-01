function [t,amps,data,aux] = read_intan_data2b(filename)
% [t,amps,data,aux] = read_intan_data
%
% Opens file selection GUI to select and then read data from an Intan 
% amplifier data file (*.int).
%
% t = time vector (in seconds)
% amps = vector listing active amplifier channels
% data = matrix of electrode-referred amplifier signals (in microvolts)
% aux = matrix of six auxiliary TTL input signals
%
% Example usage:
%  >> [t,amps,data,aux] = read_intan_data;
%  >> plot(t,data(:,1));
%
% Version 1.1, June 26, 2010
% (c) 2010, Intan Technologies, LLC
% For more information, see http://www.intantech.com
% For updates and latest version, see http://www.intantech.com/software.html
%
% 06-22-10 Added GUI file selection and optimized: Craig Patten, Plexon, Inc.

% reads second half of recording

% Open
if nargin < 1
    [file, path, filterindex] = uigetfile('*.int','Select a .int file','MultiSelect', 'off');
    filename = [path,file];
end
fid = fopen(filename, 'r');

% Read first three header bytes encoding file version
version = fread(fid, 3, 'uint8')';
if version(1) ~= 128
    error('Improper data file format.');
end
if version(2) ~= 1 || version(3) ~= 1
    warning('Data file version may not be compatible with this m-file.');
end

% Amplifier channels
amp_on = fread(fid, 64, 'uint8')';
num_amps = sum(amp_on);
amps = find(amp_on);   % create a list of amplifier channels

% Length of the data
s = dir(filename);
filesize = s.bytes;
t_count = (filesize - 67) / (num_amps * 4 + 1);
t_max = t_count / 25000;
fprintf(1, '\nData file contains %0.2f seconds of data from %d amplifier channel.\n', t_max, num_amps);

% Time
t = (t_count/2:t_count-1) / 25000;
t = t';

% Go back to the end of the header
frewind(fid);  % bof
fread(fid, 3+64, 'uint8');

% Read the first 2 tetrodes
fseek(fid,(filesize-67)/2,'cof');
data2 = fread(fid,(filesize-67)/2,'uint8=>uint8');

% extract the digital data
aux_data = data2((num_amps*4)+1:num_amps*4+1:(filesize-67)/2);

% extract individual bits
aux = [bitget(aux_data,6),bitget(aux_data,5),bitget(aux_data,4),bitget(aux_data,3),bitget(aux_data,2),bitget(aux_data,1)];
clear aux_data;

% Delete the digital data
data2((num_amps*4)+1:num_amps*4+1:(filesize-67)/2) = [];

% Convert the remaining data from bytes to single
data2 = typecast(data2,'single');

% De-mux the channels
data = zeros(t_count/2,num_amps);
for ind = 1:num_amps
    data(:,ind) = data2(ind:num_amps:length(data2));
end

% Close file
fclose(fid);