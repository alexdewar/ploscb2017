function vf_kernels2tkernels
    load('vf_kernels.mat')

    tkerns_r2 = kernel2tkernel(vf_avkernels_r2);
    tkerns_r4 = kernel2tkernel(vf_avkernels_r4);
    
    save('vf_tkernels.mat','neuroncolormap','tkerns_r2','tkerns_r4');
end
