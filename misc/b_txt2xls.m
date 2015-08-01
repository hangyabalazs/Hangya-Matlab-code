function b_txt2xls
%TXT2XLS   Converts Gabor's text files to excel files.

% File names
inpdir = 'X:\Kvantitativ_Neuroanatomia\Gabor_NO\bf\'
files = b_filelist(inpdir);
sf = length(files);

% Read text files
abc = {'A' 'B' 'C' 'D' 'E' 'F' 'G' 'H' 'I' 'J' 'K' 'L' 'M' 'N' 'O' 'P'};
for i = 1:sf
    ntxaa = files(i).name;
    cd(inpdir)
    [mtxaa mtxa] = textread(ntxaa,'%f %f');
    ntxa = {['bf_' ntxaa(1:end-4)]};
    cd ..
    eval(['xlswrite(''bf'',ntxa,''Sheet1'',''' abc{i} num2str(1) ''')']);
    eval(['xlswrite(''bf'',mtxa,''Sheet1'',''' abc{i} num2str(2) ''')']);
end