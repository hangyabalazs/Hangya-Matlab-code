function zshift_xlsprocess2
%ZSHIFT_XLSPROCESS2    Generates plots from ZSHIFT output excel file.
%   ZSHIFT_XLSPROCESS2 works on the output excel file of ZSHIFTRUN (modified
%   to contain burstiness of the theta segment with maximal z-peak for each
%   cell and information about containing parvalbumine calcium binding
%   protein). It generates various plots.
%
%   See also ZSHIFTRUN2, ZSHIFTRUN3, ZSHIFT_XLSPROCESS and ZSHIFT_XLSPROCESS3.

% Input argument check
error(nargchk(0,0,nargin))

% Read excel file
global DATAPATH
fn = [DATAPATH 'Zshift_rao\text\zshift2_theta.xls'];
mtx = xlsread(fn,'maxlen');

% Plot I - raw
figure
hold on
for i = 1:size(mtx,1)
    plot(mtx(i,3),mtx(i,4),'kx','MarkerSize',10);
end

xlabel('Z-shift')
ylabel('Z-max')

% Plot II - identified
figure
hold on
for i = 1:size(mtx,1)
    if mtx(i,10) == 1 & mtx(i,11) == 1
        h1 = plot(mtx(i,3),mtx(i,4),'mo','MarkerFaceColor','magenta');
    elseif mtx(i,10) == 1 & mtx(i,11) ~= 1
        h2 = plot(mtx(i,3),mtx(i,4),'ro','MarkerFaceColor','red');
    elseif mtx(i,10) ~= 1 & mtx(i,11) == 1
        h3 = plot(mtx(i,3),mtx(i,4),'bo','MarkerFaceColor','blue');
    elseif mtx(i,10) ~= 1 & mtx(i,11) ~= 1
        h4 = plot(mtx(i,3),mtx(i,4),'kx');
    end
end

legend([h1 h2 h3],'dual','PV','HCN')
xlabel('Z-shift')
ylabel('Z-max')

% Plot III - groups
figure
hold on
for i = 1:size(mtx,1)
    if mtx(i,12) == 1
        h1 = plot(mtx(i,3),mtx(i,4),'bo','MarkerFaceColor','blue');
    elseif mtx(i,12) == 2
        h2 = plot(mtx(i,3),mtx(i,4),'co','MarkerFaceColor','cyan');
    elseif mtx(i,12) == 3
        h3 = plot(mtx(i,3),mtx(i,4),'mo','MarkerFaceColor','magenta');
    elseif mtx(i,12) == 4
        h4 = plot(mtx(i,3),mtx(i,4),'yo','MarkerFaceColor','yellow');
    elseif mtx(i,12) == 5
        h5 = plot(mtx(i,3),mtx(i,4),'ro','MarkerFaceColor','red');
    elseif mtx(i,12) == 6
        h6 = plot(mtx(i,3),mtx(i,4),'k*','MarkerFaceColor','black');
    elseif mtx(i,12) == 7
        h7 = plot(mtx(i,3),mtx(i,4),'gp','MarkerFaceColor','green');
    else
        h8 = plot(mtx(i,3),mtx(i,4),'kx');
    end
end

legend([h1 h2 h3 h4 h5 h6 h7],'klasszikus megforduló','septum-dominált félmegforduló',...
    'hippocampus-dominált félmegforduló','hippocampus-dominált',...
    'fél-hippocampus-dominált','semmi','septum-dominált')
xlabel('Z-shift')
ylabel('Z-max')

% Plot IV - united groups
figure
hold on
for i = 1:size(mtx,1)
    if mtx(i,12) == 1 | mtx(i,12) == 2 | mtx(i,12) == 7
        h1 = plot(mtx(i,3),mtx(i,4),'ro','MarkerFaceColor','red');
    elseif mtx(i,12) == 4 | mtx(i,12) == 5
        h2 = plot(mtx(i,3),mtx(i,4),'bo','MarkerFaceColor','blue');
    else
        h8 = plot(mtx(i,3),mtx(i,4),'kx');
    end
end

legend([h1 h2],'septum vezet','hippocampus vezet')
xlabel('Z-shift')
ylabel('Z-max')

% Plot V - identified+united groups
figure
hold on
for i = 1:size(mtx,1)
    A1 = mtx(i,12) == 1 | mtx(i,12) == 2 | mtx(i,12) == 7;
    A2 = mtx(i,12) == 4 | mtx(i,12) == 5;
    B1 = mtx(i,10) == 1 & mtx(i,11) == 1;
    B2 = mtx(i,10) == 1 & mtx(i,11) ~= 1;
    B3 = mtx(i,10) ~= 1 & mtx(i,11) == 1;
    B4 = mtx(i,10) ~= 1 & mtx(i,11) ~= 1;
    if A1
        h1 = plot(mtx(i,3),mtx(i,4),'kv','MarkerSize',12);
        if B1
            set(h1,'MarkerFaceColor','magenta');
        elseif B2
            set(h1,'MarkerFaceColor','red');
        elseif B3
            set(h1,'MarkerFaceColor','blue');
        end
    elseif A2
        h2 = plot(mtx(i,3),mtx(i,4),'ks','MarkerSize',12);
        if B1
            set(h2,'MarkerFaceColor','magenta');
        elseif B2
            set(h2,'MarkerFaceColor','red');
        elseif B3
            set(h2,'MarkerFaceColor','blue');
        end
    else
        h3 = plot(mtx(i,3),mtx(i,4),'ko','MarkerSize',12);
        if B1
            set(h3,'MarkerFaceColor','magenta');
        elseif B2
            set(h3,'MarkerFaceColor','red');
        elseif B3
            set(h3,'MarkerFaceColor','blue');
        end
    end
end
set(gca,'FontSize',16)

xlabel('Z-shift')
ylabel('Z-max')