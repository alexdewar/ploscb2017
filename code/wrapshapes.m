function [nx,ny] = wrapshapes(x,y,nbins)
if nargin < 3
    nbins = 100;
end

r = max(x(:))/(2*pi);
xbin = linspace(min(x(:)),max(x(:)),nbins);
xstart = xbin(1:end-1);
xend = xbin(2:end);
[cx,cy] = deal(cell(1,size(x,2)));
for i = 1:size(x,2)
    
end