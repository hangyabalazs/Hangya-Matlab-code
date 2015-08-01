function zshift_xlsprocess
%ZSHIFT_XLSPROCESS    Generates plots from ZSHIFT output excel file.
%   ZSHIFT_XLSPROCESS works on the output excel file of ZSHIFTRUN (modified
%   to contain burstiness of the theta segment with maximal z-peak for each
%   cell and information about containing parvalbumine calcium binding
%   protein). It generates various plots.
%
%   See also ZSHIFTRUN2, ZSHIFTRUN3, ZSHIFT_XLSPROCESS2 and
%   ZSHIFT_XLSPROCESS3.

% Input argument check
error(nargchk(0,0,nargin))

% Read excel file
global DATAPATH
fn = [DATAPATH 'Zshift\Zshift3\text\zshift_theta3.xls'];
mtx = xlsread(fn);

% Plot I
figure
hold on
for i = 1:size(mtx,1)
    if mtx(i,10) == 0
        h1 = plot(mtx(i,3),mtx(i,4),'bo');
    elseif mtx(i,10) == 1
        h2 = plot(mtx(i,3),mtx(i,4),'ro');
    else
        h3 = plot(mtx(i,3),mtx(i,4),'ko');
    end
end

legend([h1 h2 h3],'non-bursting','bursting','not specified');
xlabel('Z-shift')
ylabel('Z-max')

% Plot II
figure
hold on
for i = 1:size(mtx,1)
    if mtx(i,10) == 0
        h1 = plot(mtx(i,3),mtx(i,4),'bo');
    elseif mtx(i,10) == 1
        h2 = plot(mtx(i,3),mtx(i,4),'ro');
    end
end

legend([h1 h2],'non-bursting','bursting');
xlabel('Z-shift')
ylabel('Z-max')

% Plot III
figure
hold on
for i = 1:size(mtx,1)
    if mtx(i,11) ~= 1
        h2 = plot(mtx(i,3),mtx(i,4),'kx','MarkerEdgeColor','white');
    end
end
for i = 1:size(mtx,1)
    if mtx(i,11) == 1
        h1 = plot(mtx(i,3),mtx(i,4),'ro','MarkerFaceColor','red');
    end
end
set(gca,'Color','black')

L = legend(h1,'PV-positive');
set(L,'Color','white')
xlabel('Z-shift')
ylabel('Z-max')

% Plot IV
figure
hold on
for i = 1:size(mtx,1)
    if mtx(i,11) == 1
        h1 = plot(mtx(i,3),mtx(i,4),'ro');
    end
end

legend(h1,'PV-positive');
xlabel('Z-shift')
ylabel('Z-max')

% Plot V
figure
hold on
for i = 1:size(mtx,1)
    if mtx(i,12) == 1
        h1 = plot(mtx(i,3),mtx(i,4),'ro');
    end
end

legend(h1,'PV-positive');
xlabel('Z-shift')
ylabel('Z-max')

% Plot VI
figure
hold on
for i = 1:size(mtx,1)
    if mtx(i,12) == 1
        h1 = plot(mtx(i,5)*180/pi,mtx(i,6),'ro');
    end
end

legend(h1,'PV-positive');
xlabel('Hilbert angle')
ylabel('Hilbert mvl')

% Plot VII
figure
hold on
for i = 1:size(mtx,1)
    if mtx(i,12) == 1
        h1 = plot(mtx(i,7)*180/pi,mtx(i,8),'ro');
    end
end

legend(h1,'PV-positive');
xlabel('Crosswavelet angle')
ylabel('Crosswavelet mvl')

% Plot VIII
figure
hold on
for i = 1:size(mtx,1)
    if mtx(i,12) == 1
        h1 = plot(mtx(i,5)*180/pi,mtx(i,7)*180/pi,'ro');
    end
end

legend(h1,'PV-positive');
xlabel('Hilbert angle')
ylabel('Crosswavelet angle')