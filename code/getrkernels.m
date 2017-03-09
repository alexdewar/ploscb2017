function rkernels=getrkernels(im_size,r2,r4)
if nargin < 3
    r4 = true;
end
if nargin < 2
    r2 = true;
end

load('vf_kernels.mat');
if r2
    kernels = vf_avkernels_r2;
else
    kernels = [];
end
if r4
    kernels = [kernels,vf_avkernels_r4];
end

rkernels = NaN([im_size,length(kernels)]);
for i = 1:length(kernels)
    rkernels(:,:,i) = resizekernel(kernels(i).k,im_size,thresh);
end