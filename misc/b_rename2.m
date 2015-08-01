function b_rename2

inpdir = 'D:\_analysis\matlab_data\Entry3\Theta\dtf\onesec\real\';
destdir = 'D:\_analysis\matlab_data\Entry3\Noth\dtf\onesec\real\';
mm = pwd;
cd(inpdir)
files = dir(inpdir);
sf = length(files);
for i = 1:sf
    fn = files(i).name;
    if findstr(fn,'noth')
        source = fullfile(inpdir,fn);
        copyfile(source,destdir)
        delete(fn)      
    end
end
cd(mm)