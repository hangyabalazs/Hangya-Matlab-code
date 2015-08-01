function b_rename

inpdir = 'D:\_analysis\matlab_data\Entry3\Noth\entropy\power\line\windowsize1000_overlap2\real\';
mm = pwd;
cd(inpdir)
files = dir(inpdir);
sf = length(files);
for i = 1:sf
    fn = files(i).name;
    if findstr(fn,'CLUSTER')
        cmps = strread(fn,'%s','delimiter','_');
        fn2 = [cmps{1} '_' cmps{2} '_' cmps{3} '_' cmps{4} '_' cmps{5} '_' cmps{7}];
        load(fn)
        save(fn2,'aHx','aHy','aHxy','aIxy','aIxynorm','aRelHy','aRelHx',...
            'aHxcy','aHycx','aUxy','aUyx');
        delete(fn)      
    end
end
cd(mm)