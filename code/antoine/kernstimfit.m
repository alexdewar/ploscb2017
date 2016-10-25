ktype = 'r2';
knum = 8;
thresh = 0.8;
origstim = im2double(loadim(fullfile(mfiledir,'Tpic.jpg')));

load('vf_kernels.mat',['vf_avkernels_' ktype]);
kern = eval(sprintf('vf_avkernels_%s(%d)',ktype,knum));

origstim = imresize(origstim,[NaN size(kern.k,2)]);
origstim = origstim >= thresh;
nrow = (size(kern.k,1)-size(origstim,1))/2;
ncol = size(origstim,2);
origstim = [ones(nrow,ncol); origstim; ones(nrow,ncol)];

origact = sum(sum((1-origstim).*kern.k));
origpos = sum(sum((1-origstim).*(kern.k>0)));
origneg = sum(sum((1-origstim).*(kern.k<0)));

npos = sum(sum(kern.k>0));
nneg = sum(sum(kern.k<0));
posact = 1/npos;
negact = 1/nneg;

% figure(1);clf
% subplot(1,2,1)
% showkernel(kern);
% subplot(1,2,2)
% imshow(origstim);