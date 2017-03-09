% function makegaborkernels
load('vf_kernels.mat','neuroncolormap');
thresh = 0.25;

ths = 0:pi/4:3*pi/4;
gb_kernels = cell(1,numel(ths));
for i = 1:numel(ths)
    gb = gabor_fn(1,1,0,20,ths(i));
    gb = im2kern(gb,thresh);
    
%     figure(1);clf
%     showkernel(gb);
%     keyboard
    
    gb_kernels{i} = gb;
end

save('gb_kernels.mat','gb_kernels','thresh','ths','neuroncolormap');