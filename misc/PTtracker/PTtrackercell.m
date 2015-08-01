%%

vid = videoinput('winvideo', 1);
set(vid, 'ReturnedColorSpace', 'RGB');
preview(vid)

%%

pause(10)
I = getsnapshot(vid);
H = figure;
imshow(I);

%%

imagedir = 'c:\Balazs\_personal\PTtracker\images\';
dsr = datestr(now);
dsrs = regexprep(dsr,':','_');
fname = fullfile(imagedir,['arm_' dsrs '.fig']);
% saveas(H,fname)

%%

figure(H)
k = waitforbuttonpress;
point1 = get(gca,'CurrentPoint');
point1 = point1(1,1:2);
k = waitforbuttonpress;
point2 = get(gca,'CurrentPoint');
point2 = point2(1,1:2);
k = waitforbuttonpress;
point3 = get(gca,'CurrentPoint');
point3 = point3(1,1:2);

%% angle

x = point1 - point2;
y = point3 - point2;
xy = dot(x,y);

% Arm length
lenx = sqrt(sum(x.^2));
leny = sqrt(sum(y.^2));

% Arm angle
angle = acos(xy/(lenx*leny)) * 180 / pi

%% save

datadir = 'c:\Balazs\_personal\PTtracker\data\';
fn = fullfile(datadir,['extension_' dsrs '.mat']);
% fn = fullfile(datadir,['flexion_' dsrs '.mat']);
% save(fn,'angle','dsr')

%%

line([point1(1) point2(1) point3(1)],[point1(2) point2(2) point3(2)],...
    'LineWidth',3,'Color',[0.75 0 0])
% set(gca,'Unit','normalized')
text(point2(1)+50,point2(2),['\it ' num2str(angle)],'FontSize',18,'Color',[0.75 0 0])
text
