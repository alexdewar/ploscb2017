load('kernels.mat')
kernelnum = 25;
kernelscale = 1/10;

figure(1);clf
imagesc(kernels{kernelnum});
colormap(drosophilacolormap);
colorbar
title('Kernel before thresholding')

figure(2);clf
imagesc(tkernels{kernelnum});
colormap(drosophilacolormap);
colorbar
title('Kernel after thresholding')

cm = imread('cameraman.tif');
figure(3);clf
imshow(cm);
hold on
sz = kernelscale*size(tkernels{kernelnum});
szcm = size(cm)/2;
fill(szcm(2)+sz(2)*[-1 -1 1 1 -1],szcm(1)+sz(1)*[-1 1 1 -1 -1],'r','LineStyle','none','FaceAlpha',0.5);
title('Image to process (kernel size in red)')

figure(4);clf
cmc = conv2(imresize(im2double(cm),1/kernelscale),tkernels{kernelnum});
imshow(cmc);
title('Image after convolution');