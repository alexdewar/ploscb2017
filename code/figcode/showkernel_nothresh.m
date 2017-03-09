function showkernel_nothresh(kern)
load('vf_kernels','neuroncolormap')

pos = kern > 0;
neg = kern < 0;

kern(pos) = kern(pos)/max(kern(pos));
kern(neg) = -kern(neg)/min(kern(neg));

imagesc(kern);
colormap(neuroncolormap);
