clear

load('vf_kernels_nothresh','vf_avkernels*','neuroncolormap')
kerns = [vf_avkernels_r2,vf_avkernels_r4];

figure(1);clf
for i = 1:length(kerns)
    subplot(7,6,i)
    imagesc(kerns(i).k)
    colormap(neuroncolormap)
end
