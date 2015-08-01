%%

Bcg = nan(200,480,640,3);
for k = 1:200
    Bcg(k,:) = imread('c:\Balazs\courses\2013_TENSS\tracking\mousevideo1.tif',1);
    
end
Bcg2 = median(Bcg,1);


%%

I = imread('c:\Balazs\courses\2013_TENSS\tracking\mousevideo1.tif',1);

%%

imshow(I)
figure
imagesc(I)

%%

Ic = imcrop(I,[1 1 640 410]);
figure
imagesc(Ic)

%%

Ig = rgb2gray(Ic);
figure
imshow(Ig)

%%

Ib = image2binary(Ig,125);
figure
imshow(Ib)