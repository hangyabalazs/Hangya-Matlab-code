%%

load('C:\Balazs\_analysis\NB\regression\All\regression_results.mat')

%%

% load('C:\Balazs\_analysis\NB\regression\All\new2\regression_results.mat')
load('C:\Balazs\_analysis\NB\regression\All\impexclude\regression_results.mat')
NumCells = length(p);
[p2 R2] = deal(nan(1,NumCells));
for k = 1:NumCells
    if ~isempty(p(k).iti05)
        p2(k) = p(k).iti05;
        R2(k) = R(k).iti05;
    end
end
p = p2;
R = R2;

%%

% area = getvalue('Area1');
% pNBinx = cellfun(@(s)ismember(s,{'GP','GP/SI','SI','IC','RT/IC','EP','EA','EAC'}),area);
% ID_PC = getvalue('ID_PC');
% Lr_PC = getvalue('Lr_PC');
% validity = getvalue('validity');

NB = selectcell(['"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})']);
[~, NBinx] = intersect(I,NB);

VB = selectcell(['"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("Area1",{''VPM'',''VPL'',''VB''})']);
[~, VBinx] = intersect(I,VB);

CPu = selectcell(['"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("Area1",{''CPu'',''Cpu''})']);
[~, CPuinx] = intersect(I,CPu);

RT = selectcell(['"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("Area1",{''RT''})']);
[~, RTinx] = intersect(I,RT);

ACx = selectcell(['"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("Area1",{''ACx''})']);
[~, ACxinx] = intersect(I,ACx);

%%

ChAT = selectcell(['"ChAT+"==1&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})']);  % identified ChAT+ cells
ChAT = [ChAT 'n045_121217x_4.6'];   % clustered based on light-evoked spikes
[~, ChATinx] = intersect(I,ChAT);

pChAT = selectcell(['"pChAT+"==1&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})']);  % putative ChAT+ cells
[~, pChATinx] = intersect(I,pChAT);

PV = selectcell(['"PV+"==1&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})']);  % identified PV+ cells
[~, PVinx] = intersect(I,PV);


%%

nninx = find(~isnan(p));
NBbinx = intersect(NBinx,nninx);
NBsinx = intersect(NBinx,find(p<0.01));
CORR = I(NBsinx);

%%

sum(p(NBinx)<0.01)/sum(~isnan(p(NBinx)))
sum(p(RTinx)<0.01)/sum(~isnan(p(RTinx)))
sum(p(VBinx)<0.01)/sum(~isnan(p(VBinx)))
sum(p(ACxinx)<0.01)/sum(~isnan(p(ACxinx)))
sum(p(CPuinx)<0.01)/sum(~isnan(p(CPuinx)))

%%

sum(p(ChATinx)<0.01)/sum(~isnan(p(ChATinx)))
sum(p(pChATinx)<0.01)/sum(~isnan(p(pChATinx)))
sum(p(PVinx)<0.01)/sum(~isnan(p(PVinx)))


%%

PSTH = getvalue('Hit_psth');
PSTH = cell2mat(PSTH(NBsinx));

%%

PSTHn = nan(size(PSTH));
for k = 1:size(PSTH,1)
    PSTHn(k,:) = zscore(PSTH(k,:));
end

%%

figure;plot(mean(PSTHn))
figure;imagesc(PSTHn)

%% call 'nbstimraster'

nbstimraster(I(NBsinx),true)