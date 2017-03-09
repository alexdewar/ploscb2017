function showkernel(kern)
exc = [1 0 0];
inc = [0 0 1];

if isstruct(kern)
    kern = kern.k;
end
imagesc(sign(kern));
colormap([inc;1 1 1;exc]);