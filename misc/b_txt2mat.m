function b_txt2mat
%TXT2MAT   Converts txt 3-channel data file to mat file.
%   TXT2MAT reads 3-channel data file (unit-unit registration) from text
%   file and saves mat file in the same directopry. Specify input directory
%   via editing the program code!
%
%   See also TEXTREAD.

% Input argumnet check
error(nargchk(0,0,nargin))

% Input directory
mm = pwd;
inpdir = ['f:\raw_data\unit_unit\temp\'];
cd(inpdir)

% Read and write
files = dir(inpdir);
sf = length(files)
for i = 1:sf;
    [fpath fname fext] = fileparts(files(i).name);
    if isequal(fext,'.txt')
        [A B C D] = textread([fname '.txt'],'%f %f %f %f','headerlines',2);
        data = [B C D];
        save(fname,'data')
    end
end
cd(mm)    