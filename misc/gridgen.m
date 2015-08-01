function gridgen
%GRIDGEN   Grid cell rate map generation.
%   Based on the idea of Alexei V. Samsonovich

n = 200;
m = 100;
mh = 70;

% Generate random local view input for n=200 channels
x = randn(m,m,n);

% Smoothen input channels as functions of maze location
h = exp(-((cumsum(ones(mh),1)-mh/2).^2+(cumsum(ones(mh),2)-mh/2).^2)./2./mh);
h = h(1:end-1,1:end-1);
for k = 1:n
    y(:,:,k) = filter2(h,x(:,:,k));
end

% Convert the vector of input channels to principal components
y2 = reshape(y,m*m,n);
[pc,y2] = princomp(y2);
y2 = reshape(exp(y2),m,m,n);

% Plot the result
figure
pcolor(y2(:,:,n));
axis equal
% colormap cool
shading flat