function b_entropytrf
%ENTROPYTRF    Rename files.

% Import
global DATAPATH
where1 = [DATAPATH,'Entropy\entropy_20__5000\temp2\'];    %Here are the source files
files = dir(where1);
files = files(3:end);
sf = length(files);
mm = pwd;
cd([DATAPATH,'Entropy\entropy_20__5000\temp2\']);  %Here are the destination files

wb = waitbar(0,'Running ENTROPYTRF...','Position',[360 250 275 50]);    %Progress indicator
for o = 1:sf
    if ~files(o).isdir
        fname = files(o).name;
        ffnm = [where1 fname];
        load(ffnm);
        cl = fname(1:end-4);
        str = [cl '_ENTROPY'];
        save(str,'Rhxabs','Rhyabs','Rhxyabs','Rixyabs','Rixynormabs','Rrelshanabs','Rhxcyabs',...
            'Rhycxabs','Ruxyabs','Ruyxabs');
    end
    
    waitbar(o/sf)   %Progress indicator
end
close(wb);   %Close progress indicator