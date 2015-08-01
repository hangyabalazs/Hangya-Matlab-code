%% load

load('C:\MATLAB_R2010a\work\Adam\Psth_Matrix_HomeZoneOut1.mat')
non_tagged_MatrixPsth(29,:) = [];   % all NaNs

%% load - new data

load('C:\My Dropbox\KepecsLab\_Duda\for balazs\HomeOutRasterMatrix.mat')
non_tagged_MatrixPsth = Non_Tagged_Matrix;
pv_MatrixPsth = PV_Matrix;
som_MatrixPsth = SOM_Matrix;

%% test if overlapping

all = [pv_MatrixPsth; som_MatrixPsth; non_tagged_MatrixPsth];
isequal(all(14,:),all(171,:))
figure
plot(all(2,:))
hold on
plot(all(2,:),'r')

%% plot

figure;plot(pv_MatrixPsth','r')
hold on;plot(som_MatrixPsth','b')

%% get rid of edge effect

pv_MatrixPsth = pv_MatrixPsth(:,10:end-10);
som_MatrixPsth = som_MatrixPsth(:,10:end-10);
% non_tagged_Inh_MatrixPsth = non_tagged_Inh_MatrixPsth(:,10:end-10);
non_tagged_MatrixPsth = non_tagged_MatrixPsth(:,10:end-10);

%% smooth

sizepv = size(pv_MatrixPsth);
spv = zeros(sizepv(1),sizepv(2));
for k = 1:sizepv(1)
    spv(k,:) = smooth(pv_MatrixPsth(k,:),'linear',101);
end

sizesom = size(som_MatrixPsth);
ssom = zeros(sizesom(1),sizesom(2));
for k = 1:sizesom(1)
    ssom(k,:) = smooth(som_MatrixPsth(k,:),'linear',101);
end

% sizeinh = size(non_tagged_Inh_MatrixPsth);
% sinh = zeros(sizeinh(1),sizeinh(2));
% for k = 1:sizeinh(1)
%     sinh(k,:) = smooth(non_tagged_Inh_MatrixPsth(k,:),'linear',101);
% end

sizent = size(non_tagged_MatrixPsth);
snt = zeros(sizent(1),sizent(2));
for k = 1:sizent(1)
    snt(k,:) = smooth(non_tagged_MatrixPsth(k,:),'linear',101);
end

figure;plot(spv','r')
hold on;plot(ssom','b')
plot(median(spv),'r','LineWidth',3)
plot(median(ssom),'b','LineWidth',3)

%% plot non-tagged

figure
plot(sinh','k')

figure
plot(snt','k')

%% distance from median

medpv = median(spv);
spvd = zeros(1,sizepv(1));
for k = 1:sizepv(1)
    spvd(k) = sum((spv(k,:)-medpv).^2);
end

medsom = median(ssom);
ssomd = zeros(1,sizesom(1));
for k = 1:sizesom(1)
    ssomd(k) = sum((ssom(k,:)-medpv).^2);
end

%% PCA

pv_MatrixPsth2 = spv(:,100:400);
som_MatrixPsth2 = ssom(:,100:400);

[coeff scores] = princomp([pv_MatrixPsth2; som_MatrixPsth2]);
PC1 = scores(:,1);
PC2 = scores(:,2);
PC3 = scores(:,3);
PC4 = scores(:,4);

figure
plot(PC3(1:sizepv(1)),PC2(1:sizepv(1)),'b.','MarkerSize',20)
hold on
plot(PC3(sizepv(1)+1:end),PC2(sizepv(1)+1:end),'r.','MarkerSize',20)

%% restrict in time

spv2 = spv(:,150:300);
ssom2 = ssom(:,150:300);
sinh2 = sinh(:,150:300);
snt2 = snt(:,150:300);

%% restrict in time

spv2 = spv(:,100:400);
ssom2 = ssom(:,100:400);
% sinh2 = sinh(:,100:400);
snt2 = snt(:,100:400);

%% restrict in time - new data

spv2 = spv(:,100:500);
ssom2 = ssom(:,100:500);
snt2 = snt(:,100:500);

%% PCA - all included

[coeff scores] = princomp([spv2; ssom2; snt2]);
PC1 = scores(:,1);
PC2 = scores(:,2);
PC3 = scores(:,3);
PC4 = scores(:,4);

figure
plot(PC3(sizepv(1)+sizesom(1)+1:end),PC2(sizepv(1)+sizesom(1)+1:end),'.','MarkerSize',20,'Color',[0.7 0.7 0.7])
hold on
plot(PC3(1:sizepv(1)),PC2(1:sizepv(1)),'b.','MarkerSize',20)
plot(PC3(sizepv(1)+1:sizepv(1)+sizesom(1)),PC2(sizepv(1)+1:sizepv(1)+sizesom(1)),'r.','MarkerSize',20)

%% Voronoi diagram

x = [PC3(1:sizepv(1)+sizesom(1)) PC2(1:sizepv(1)+sizesom(1))];
% x = [PC3 PC2];
[v,c] = voronoin(x); 
for pg = 1:length(c)
    if all(c{pg}~=1)   % If at least one of the indices is 1, then it is an open region and we can't patch that.
        ip = inpolygon(x(:,1),x(:,2),v(c{pg},1),v(c{pg},2));
        if ismember(find(ip),1:sizepv(1))
            edgeclr = [0 0 1];
            faceclr = [0.25 0.25 0.75];
        elseif ismember(find(ip),sizepv(1)+1:sizepv(1)+sizesom(1))
            edgeclr = [1 0 0];
            faceclr = [0.75 0.25 0.25];
        end
        patch(v(c{pg},1),v(c{pg},2),edgeclr,'FaceColor',faceclr,'FaceAlpha',0.5); % use color i.
    end
end

%% Voronoi diagram on standardized PCA components

PC3 = standardize(PC3);
PC2 = standardize(PC2);
figure
plot(PC3(sizepv(1)+sizesom(1)+1:end),PC2(sizepv(1)+sizesom(1)+1:end),'.','MarkerSize',20,'Color',[0.7 0.7 0.7])
hold on
x = [[PC3(1:sizepv(1)+sizesom(1)); -100; -100; 100; 100]...
    [PC2(1:sizepv(1)+sizesom(1)); -100; 100; -100; 100]];
% x = [PC3 PC2];
voronoi(x(:,1),x(:,2),'k')
plot(PC3(1:sizepv(1)),PC2(1:sizepv(1)),'b.','MarkerSize',20)
plot(PC3(sizepv(1)+1:sizepv(1)+sizesom(1)),PC2(sizepv(1)+1:sizepv(1)+sizesom(1)),'r.','MarkerSize',20)

[v,c] = voronoin(x); 
for pg = 1:length(c)
    if all(c{pg}~=1)   % If at least one of the indices is 1, then it is an open region and we can't patch that.
        ip = inpolygon(x(:,1),x(:,2),v(c{pg},1),v(c{pg},2));
        if ismember(find(ip),1:sizepv(1))
            edgeclr = [0 0 1];
            faceclr = [0.25 0.25 0.75];
        elseif ismember(find(ip),sizepv(1)+1:sizepv(1)+sizesom(1))
            edgeclr = [1 0 0];
            faceclr = [0.75 0.25 0.25];
        end
        patch(v(c{pg},1),v(c{pg},2),edgeclr,'FaceColor',faceclr,'FaceAlpha',0.5); % use color i.
    end
end
axis([-4 4 -4 4])

%% clustering

fcns = [spv2; ssom2; snt2];
gr1 = 1:sizepv(1);
gr2 = sizepv(1)+1:sizepv(1)+sizesom(1);
somclust(fcns,gr1,gr2)

%% rising edge slope - peak loc. plane

nm_cells = size(fcns,1);
slope = nan(1,nm_cells);
cd('c:\Balazs\_analysis\SOM\slope\')
for k = 1:nm_cells
    [H,slope(k)] = somslope(fcns(k,:));
    saveas(H,['slope_' num2str(k)])
    close(H)
end

peaks = nan(1,nm_cells);
for k = 1:nm_cells
    pk = find(fcns(k,:)==max(fcns(k,:)));
    peaks(k) = pk(1);
end

figure
plot(slope(sizepv(1)+sizesom(1)+1:end),peaks(sizepv(1)+sizesom(1)+1:end),'.','MarkerSize',20,'Color',[0.7 0.7 0.7])
hold on
plot(slope(1:sizepv(1)),peaks(1:sizepv(1)),'b.','MarkerSize',20)

%% rising edge slope - peak loc. plane; with somslope2

nm_cells = size(fcns,1);
slope = nan(1,nm_cells);
cd('c:\Balazs\_analysis\SOM\slope2\')
for k = 1:nm_cells
    [H,slope(k)] = somslope2(fcns(k,:),80,130);
    saveas(H,['slope_' num2str(k)])
    close(H)
end

peaks = nan(1,nm_cells);
for k = 1:nm_cells
    pk = find(fcns(k,:)==max(fcns(k,1:300)));
    peaks(k) = pk(1);
end

figure
plot(slope(sizepv(1)+sizesom(1)+1:end),peaks(sizepv(1)+sizesom(1)+1:end),'.','MarkerSize',20,'Color',[0.7 0.7 0.7])
hold on
plot(slope(1:sizepv(1)),peaks(1:sizepv(1)),'b.','MarkerSize',20)
plot(slope(sizepv(1)+1:sizepv(1)+sizesom(1)),peaks(sizepv(1)+1:sizepv(1)+sizesom(1)),'r.','MarkerSize',20)