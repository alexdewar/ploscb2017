function kern=im2kern(im,thresh)
kern = sign(im).*(abs(im)>=thresh);
kern = sign(kern);
kern(kern==1) = 1./sum(kern(:)==1);
kern(kern==-1) = -1./sum(kern(:)==-1);