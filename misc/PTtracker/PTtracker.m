function PTtracker(issave)
%PTTRACKER   Track physiscal therapy progress.
%   PTTRACKER measures and stores the angle of the joint motion.

if nargin < 1
    issave = true;   % save image and angle data
end
main('extension',issave)   % extension
main('flexion',issave)   % flexion

% -------------------------------------------------------------------------
function main(motion,issave)

% Start video
vid = videoinput('winvideo', 1);   % video object
set(vid, 'ReturnedColorSpace', 'RGB');
PV = preview(vid);   % launch video

% Take a snapshot
pause(10)
I = getsnapshot(vid);   % snapshot
delete(PV)   % stop video
H = figure;
imshow(I);   % show snapshot

% Add anchor points
figure(H)
k = 1;
while k ~= 0
    k = waitforbuttonpress;
end
point1 = get(gca,'CurrentPoint');   % shoulder
point1 = point1(1,1:2);

k = 1;
while k ~= 0
    k = waitforbuttonpress;
end
point2 = get(gca,'CurrentPoint');   % elbow
point2 = point2(1,1:2);

k = 1;
while k ~= 0
    k = waitforbuttonpress;
end
point3 = get(gca,'CurrentPoint');   % wrist
point3 = point3(1,1:2);

% Angle
x = point1 - point2;
y = point3 - point2;
xy = dot(x,y);

lenx = sqrt(sum(x.^2));   % arm length
leny = sqrt(sum(y.^2));

angle = acos(xy/(lenx*leny)) * 180 / pi    % arm angle

% Put the angle in the image
line([point1(1) point2(1) point3(1)],[point1(2) point2(2) point3(2)],...
    'LineWidth',3,'Color',[0.75 0 0])
% set(gca,'Unit','normalized')
text(point2(1)+50,point2(2),['\it ' num2str(angle)],'FontSize',18,'Color',[0.75 0 0])

% Store image
imagedir = 'c:\Balazs\_personal\PTtracker\images\';
dsr = datestr(now);
dsrs = regexprep(dsr,':','_');
fname = fullfile(imagedir,['arm_' dsrs '.fig']);
if issave
    saveas(H,fname)
end

% Save angle data
datadir = 'c:\Balazs\_personal\PTtracker\data\';
fn = fullfile(datadir,[motion '_' dsrs '.mat']);
% fn = fullfile(datadir,['flexion_' dsrs '.mat']);
if issave
    save(fn,'angle','dsr')
    disp('Angle saved.')
else
    disp('Angle not saved.')
end