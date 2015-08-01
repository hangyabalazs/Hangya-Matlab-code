function cbsloperun
%CBSLOPERUN   Runs CBSLOPE on a sequence of files.
%   CBSLOPERUN calculates slopes and amplitudes for stimulus response on a
%   sequence of files; see CBSLOPE for detailed help. Edit program code to
%   modify input and output directories.
%
%   See also CBSLOPE and LINEFIT.

% Directories
inpdir = 'X:\In_Vivo\raw_data\hippo_nape_matfiles\processing';
tdir = 'D:\_analysis\hippo_nape\result_figs\';
mm = pwd;

% File list
files = filelist(inpdir);
sf = length(files);

% Call CBSLOPE
slope = cell(1,sf);
amp = cell(1,sf);
seg = zeros(1,sf);
st = zeros(1,sf);
en = zeros(1,sf);
for o = 1:sf
    cd(inpdir)
    fname = files{o};
    load(fname)
    [slope{o} amp{o}] = cbslope(data);
    cmps = strread(fname,'%s','delimiter','_');
    seg(o) = str2num(cmps{7});
    st(o) = str2num(cmps{8});
end

% Generate output
Slope = [];
Amp = [];
while ~isempty(seg)
    seginx = find(seg==min(seg));
    stnew = st(seginx);
    [sst ixst] = sort(stnew);
    for k = 1:length(ixst)
        for m = 1:length(slope{seginx(ixst(k))})
            Slope = [Slope slope{seginx(ixst(k))}(m)];
            Amp = [Amp amp{seginx(ixst(k))}(m)];
        end
    end
    seg(seginx) = [];
    st(seginx) = [];
    slope(seginx) = [];
    Slope = [Slope NaN];
    amp(seginx) = [];
    Amp = [Amp NaN];
end

% Plot ans save result
figure(1)
plot(Slope,' k.')
naninx = find(isnan(Slope));
y_lim = ylim;
line([naninx; naninx],[repmat(y_lim(1),1,length(naninx)); repmat(y_lim(2),1,length(naninx))],...
    'Color','red')

saveas(gcf,[tdir fname(1:3) '_slope.tif'],'tiffn');

figure(2)
plot(Amp,' k.')
naninx = find(isnan(Amp));
y_lim = ylim;
line([naninx; naninx],[repmat(y_lim(1),1,length(naninx)); repmat(y_lim(2),1,length(naninx))],...
    'Color','red')

saveas(gcf,[tdir fname(1:3) '_amp.tif'],'tiffn');

cd(mm)



% -------------------------------------------------------------------------
function files2 = filelist(inpdir)

files = dir(inpdir);
files = files(3:end);
files2 = {};
for i = 1:length(files)
    if ~files(i).isdir
        files2{end+1} = files(i).name;
    end
end