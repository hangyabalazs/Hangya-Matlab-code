function ainverteeg
%AINVERTEEG   Invert EEG of all files in a directory.
%   AINVERTEEG loads all files in a specified directory, inverts EEG and
%   overwrites the loaded file with the new data.

inpdir = 'F:\balazs\_analysis\Andi\Disc\ketxyl\control\AUJ31\'
cd(inpdir)
files = dir(pwd);
for k = 1:length(files)
    if files(k).isdir
        continue
    end
    ff = files(k).name;
    load(ff)
    eeg = data(:,2);
    eeg = -eeg;
    data = [data(:,1) eeg];
    save(ff,'data')
    clear eeg data ff
end