%% delete CSC files

fs = filesep;

ap = 'F:\NB_cellbase\n026\';
dra = dir(ap);

for k = 3:length(dra)
    sp = [ap dra(k).name];
    fn1 = [sp fs 'CSC1.ncs'];
    fn2 = [sp fs 'CSC2.ncs'];
    fn3 = [sp fs 'CSC3.ncs'];
    fn4 = [sp fs 'CSC4.ncs'];
    fn6 = [sp fs 'CSC6.ncs'];
    fn7 = [sp fs 'CSC7.ncs'];
    delete(fn1,fn2,fn3,fn4,fn6,fn7)
end

%% copy ntt files

fs = filesep;
xroot = 'x:\Data\balazs\n021\';
froot = 'f:\NB_cellbase\n021\';

dra = dir(xroot);

for k = 3:length(dra)

    sp = dra(k).name;
    fn1 = [xroot sp fs 'TT1.ntt'];
    fn2 = [xroot sp fs 'TT2.ntt'];
    fn3 = [xroot sp fs 'TT3.ntt'];
    fn4 = [xroot sp fs 'TT4.ntt'];
    fn5 = [xroot sp fs 'TT5.ntt'];
    fn6 = [xroot sp fs 'TT6.ntt'];
    fn7 = [xroot sp fs 'TT7.ntt'];
    fn8 = [xroot sp fs 'TT8.ntt'];
    dest = [froot sp];
    
    copyfile(fn1,dest)
    copyfile(fn2,dest)
    copyfile(fn3,dest)
    copyfile(fn4,dest)
    copyfile(fn5,dest)
    copyfile(fn6,dest)
    copyfile(fn7,dest)
    copyfile(fn8,dest)
end