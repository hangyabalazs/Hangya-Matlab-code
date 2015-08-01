function czspikewidth
%CZSPIKEWIDTH   Spike width calculation.
%   CZSPIKEWIDTH calculates action potential width at 1/e of the peak
%   amplitute, mean firing rate, minimal and maximal interspike interval. 
%   Results are saved in an Excel file.
%
%   See also CZPHASE2.

% Directories
global DATAPATH
inpdir_unit = [DATAPATH 'Czurko\discriminated2\'];
inputxls = [DATAPATH 'Czurko\spike_width\cluster_info.xls'];
outputxls = [DATAPATH 'Czurko\spike_width\spike_width.xls'];

% Load
headerrows = 1;
intpyr = 'PYR';
[mtx ntx atx] = xlsread(inputxls,intpyr);   % load zshift
ntx(1:headerrows,:) = [];
atx(1:headerrows,:) = [];
xlsout = atx;
fns = atx(:,1);

for o = 1:length(fns)
    cell_code = atx{o,2};
    
    % Import
    dr = [inpdir_unit atx{o,3} '\' atx{o,1} '\' atx{o,4} '\'];
    fn = [dr 'clusters.mat'];
    load(fn)        % load unit

    % Cluster indeces
    inx = find(CellNumbers==cell_code);
    
    % Find the shank with maximal amplitude
    pspikes = cell(1,4);
    pmx = cell(1,4);
    pspikes{1} = squeeze(DataPoints(:,1,inx));
    pmx{1} = max(pspikes{1},[],1);
    mmx1 = mean(pmx{1});
    pspikes{2} = squeeze(DataPoints(:,2,inx));
    pmx{2} = max(pspikes{2},[],1);
    mmx2 = mean(pmx{2});
    pspikes{3} = squeeze(DataPoints(:,3,inx));
    pmx{3} = max(pspikes{3},[],1);
    mmx3 = mean(pmx{3});
    pspikes{4} = squeeze(DataPoints(:,4,inx));
    pmx{4} = max(pspikes{4},[],1);
    mmx4 = mean(pmx{4});
    mmxs = [mmx1 mmx2 mmx3 mmx4];
    minx = find(mmxs==max(mmxs));
    spikes = pspikes{minx};
    mx = pmx{minx};
    
    % Spike width
    mxe = mx / exp(1);
    spno = length(mxe);
    srwf = 30303;      % sampling rate of waveforms
    SpikeWidth = zeros(1,spno);
    for k = 1:spno
        spk = spikes(:,k);  % current spike
        mxl = find(spk==max(spk));
        mxl = mxl(1);
        spk_h1 = spk(1:mxl);    % spike first half (before peak)
        spk_h2 = spk(mxl+1:end);    % spike second half (after peak)
        fp1 = find(spk_h1<mxe(k),1,'last');
        fp2 = find(spk_h2<mxe(k),1,'first');
        if isempty(fp1) || isempty(fp2)
            SpikeWidth(k) = NaN;
            continue
        end
        tmxe1 = fp1 + (mxe(k) - spk(fp1)) / (spk(fp1+1) - spk(fp1));   % interpolate exact 1/e crossing
        tmxe2 = mxl + fp2 - (mxe(k) - spk(mxl+fp2)) / (spk(mxl+fp2-1) - spk(mxl+fp2));   % interpolate exact 1/e crossing
        SpikeWidth(k) = (tmxe2 - tmxe1) / srwf * 10^6;  % spike width in microsec
    end
    xlsout{o,5} = nanmean(SpikeWidth);
    xlsout{o,6} = nanstd(SpikeWidth);
    xlsout{o,7} = nanstd(SpikeWidth) / sqrt(length(SpikeWidth(~isnan(SpikeWidth))));
%     if nanmean(SpikeWidth) < 200
%         keyboard
%     end
    
    % Firing rate
    len = (TimeStamps(end) - TimeStamps(1)) / 1000000;
    FiringRate = spno / len;
    xlsout{o,9} = nanmean(FiringRate);
    
    % Complex-spike bursts
    vdisc = TimeStamps(inx) / 1000000;
    isi = diff(vdisc);
    MinISI = min(isi) * 1000;  % minimal ISI in ms
    MaxISI = max(isi);  % maximal ISI in s
    xlsout{o,10} = MinISI;
    xlsout{o,11} = MaxISI;
end

% Write output file
xlswrite(outputxls,xlsout,intpyr)