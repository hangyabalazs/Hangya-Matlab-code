function aclustout

% Directories
inpdir1 = '';   % cluster analysis mat files
inpdir2 = '';   % clustercut ('dec')
inpdir3 = '';   % aftercut ('ac')
inpdir4 = '';   % discriminated unit
resdir = '';
mm = pwd;

% Filelist
[files1 files_short1] = filelist(inpdir1);
[files2 files_short2] = filelist(inpdir2);
sf = length(files_short1);

% Corrigate cluster analysis results
for o = 1:sf
    burst = Burst{dec};
    after = After{dec};
    ast = sort(after);
    for k = 1:ac
        inx = find(after==ast(k));
        fnd = find(burst<inx,1,'last');
        burst(fnd) = burst(fnd+1);
    end
    
% Intraburst intervals
    ssi = vdisc;    %ssi is for "single spikes included": it will contain the first spikes of bursts
                    %and the single spikes as well
    intraburstiv = [];
    burstnum = size(burst,2);
    intraburstnum = zeros(1,sb2);
    isimtx = {};
    for j = 1:burstnum    %computing intraburstiv and allfirstspike
        b = vdisc(burst(1,j):burst(2,j));
        intraburstiv = [intraburstiv diff(b)];                 
        ssi(burst(1,j)+1:burst(2,j)) = 0;
        intraburstnum(j) = length(b);   %intraburst spike number
        for k = 1:length(b) - 1
            isimtx{length(b),k} = [isimtx(length(b),k) vdisc(burst(1,j)+k)-vdisc(burst(1,j)+k-1)];
        end
    end
    burstlength = (vdisc(burst(2,:)) - vdisc(burst(1,:))) / 20000;
    intraburstfreq = (intraburstnum - 1) ./ burstlength;
    
% Interburst intervals
    fsp2 = fsp(2:end);  %computing interburstiv
    lsp2 = lsp(1:end-1);
    interburstiv = fsp2 - lsp2;
    
% Inter-1st-sike intervals
    interfirstspike = diff(fsp);
    if ~isempty(interfirstspike)
        firstspikefreq = 20000 / mean(interfirstspike);
        firstspikefreq = (burstnum - 1)  / (vdisc(burst(2,end)) - vdisc(burst(1,1)));
    else
        firstspikefreq = NaN;
    end
    
% Beforefirst intervals
    if ~isequal(burst(1,1),1)
        beforefirstiv = ivs(burst(1,:)-1);
    else
        beforefirstiv = ivs(burst(1,2:end)-1);
    end
    Silence(dec) = min(beforefirstiv) / 20000;
    gap = (min(beforefirstiv) - max(intraburstiv)) / 20000;
    
% Afterlast intervals
    if ~isequal(burst(end,end),length(vdisc))
        afterlastiv = ivs(burst(2,:));
    else
        afterlastiv = ivs(burst(2,1:end-1));
    end
    after = afterlastiv / 20000;
    
% Burst parameters
    Burstiness = (length(intraburstiv) + burstnum) / length(vdisc);
    IntraBurstFrequency.mean = mean(intraburstfreq);
    IntraBurstFrequency.sd = std(intraburstfreq);
    IntraBurstFrequency.all = intraburstfreq;
    IntraBurstSpikeNumber.mean = mean(intraburstnum);
    IntraBurstSpikeNumber.sd = std(intraburstnum);
    IntraBurstSpikeNumber.all = intraburstnum;
    BurstLength.mean = mean(burstlength);
    BurstLength.sd = std(burstlength);
    BurstLength.all = burstlength;
    BurstFrequency = firstspikefreq;
    for x = 1:size(isimtx,1)
        for y = 1:size(isimtx,2)
            IsiMatrix.mean(x,y) = mean(isimtx{x,y});
            IsiMatrix.sd(x,y) = std(isimtx{x,y});
            IsiMatrix.num(x,y) = length(isimtx{x,y});
        end
    end

% Save excel file    
    
end