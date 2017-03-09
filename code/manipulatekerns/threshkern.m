function kern2=threshkern(kern,thresh)
kern2 = zeros(size(kern));
abvthresh = abs(kern)>=thresh;
kern2(abvthresh) = sign(kern(abvthresh));
kern2(kern2==1) = 1/sum(kern2(:)==1);
kern2(kern2==-1) = -1/sum(kern2(:)==-1);