function rapheconfocal
%RAPHECONFOCAL    Pixel counts for confocal images.
%   RAPHECONFOCAL loads confocal images and calculates the number of
%   green/red pixels in each column. Only pixels exceeding an intensity
%   threshold (in the function of the intensity distribution) are included.
%   Pixels where the intensity of both colors reaches threshold are also
%   counted.
%
%   The program is designed to evaluate the layer specificity of
%   channel-rhodopsin labeled raphe inputs and labeled serotonergic inputs
%   in the hippocampus. Loaded images contain a panel traversing the entire
%   CA1. Horizontal axis of the image is perpendicular to the surface of
%   the hippocampus.
%
%   See also RAPHECONFOCAL_RUN.

% Load images
[Igreen,cmap] = imread('X:\In_Vivo\balazs\_analysis\Raphe\zs7\zs7g_900_s02_01_l_C001','tif');
[Ired,cmap] = imread('X:\In_Vivo\balazs\_analysis\Raphe\zs7\zs7g_900_s02_01_l_C002','tif');
Igreen = double(Igreen);
Ired = double(Ired);
Igreen = Igreen / max(Igreen(:));
Ired = Ired / max(Ired(:));

% Overlay image
Iov = Igreen + Ired;

% Show images
figure
imagesc(Ired)
figure
imagesc(Igreen)
figure
imagesc(Iov)

% Remove background and count pixels
[Ig Igo] = processcolor(Igreen,2);    % remove background
[Ir Iro] = processcolor(Ired,1);
oIgo = Igo;     % overlay
oIgo(Iro==0) = 0;
oIro = Iro;
oIro(Igo==0) = 0;
Io = zeros(size(Ired));
Io(:,:,1) = oIro;
Io(:,:,2) = oIgo;
figure
imagesc(Io)

figure;     % count pixels
sIg = sum(Igo);
plot(sIg)
xlim([0 length(sIg)])
figure;
sIr = sum(Iro);
plot(sIr)
xlim([0 length(sIr)])
figure
soIr = sum(oIro);
plot(soIr,'y')
xlim([0 length(soIr)])



% ------------------------------------------------------------------------
function [I Ig] = processcolor(I,ch)

% Remove background
Ig = I(:,:,ch);
% [nm xout] = hist(double(Ig(:)),[1:256]);
% figure
% bar(xout,nm)
% cnm = cumsum(nm);
% figure
% plot(cnm)

sig = sort(Ig(:));
lig = length(Ig(:));
levs = [0.95 0.99 0.995 0.999];  % significance levels
lle = length(levs);
sl = zeros(lle,2);
for k = 1:lle
    sl(k,1) = levs(k);
    inx = ceil(lig*levs(k));
    sl(k,2) = sig(inx);     % critical values
end

Ig(Ig<sl(3,2)) = 0;
Ig(Ig>=sl(3,2)) = 1;
I(:,:,ch) = Ig;
figure
imagesc(I)