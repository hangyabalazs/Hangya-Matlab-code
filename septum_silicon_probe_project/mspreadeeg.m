function mspreadeeg
%MSPREADEEG   Import EEG data.
%   MSPREADEEG imports EEG data of the MSHCsp project and exports pyramidal
%   layer EEG to .mat files.
%
%   See also MSPACORR.

% Input argument check
error(nargchk(0,0,nargin))

% Directories
global DATAPATH
global DATADIR
inpdir = [DATADIR 'MSHCsp\Viktor5\'];
resdir = [DATADIR 'MSHCsp\Viktor5\EEG\'];
mm = pwd;

% Filelist
% rid = '1234';
% flist = {'201003171' '201003172' '201003173' '201003174'};
rid = '1237';
flist = {'201003291' '201003292' '201003293' '201003294'};
sf = length(flist);

% Main
sr = 20000;
% Inx = [33:65];
Inx = [60];
for o = 1:sf
    ff = [inpdir flist{o} '.par'];
    pars = LoadPar(ff);
    pchs = cellfun(@(x) x',pars.ElecGp,'UniformOutput',0);
    chs = cell2mat(pchs);
    
    fname = [inpdir flist{o} '.dat'];
    numch = pars.nChannels;
    buffersize = 4096*10;
    fileinfo = dir(fname);  % get file size, and calculate the number of samples per channel
    flen = ceil(fileinfo(1).bytes/2/numch);
    numelmax = flen;
    numelmin = flen - 1000000;
    numelmin = 0;
    datafile = fopen(fname,'r');
    eeg = cell(1,numch);
    for k = Inx
        eeg{k} = zeros(numelmax-numelmin+1,1);
    end
    numel = 0;
    fseek(datafile,numelmin*numch*2,'bof');
    while numelmin + numel < numelmax
        [data,count] = fread(datafile,[numch,buffersize],'int16');
        numelm = count / numch;
        for k = Inx
            chselect = chs(k) + 1;
            eeg{k}(numel+1:numel+numelm) = data(chselect,:);
        end
        numel = numel + numelm;
    end
    fclose(datafile);
    eeg = eeg{60};
    fn = [resdir rid '_' flist{o} '_eeg.mat'];
    save(fn,'eeg')
end
% figure
% for k = Inx
%     clr = rand(1,3);
%     plot(eeg{k}-mean(eeg{k})+k*2000,'Color',clr)
%     hold on
% end