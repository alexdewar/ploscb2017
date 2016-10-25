fov = 270;

im = im2double(rgb2gray(imread('./dispatchrecap/antoinestim/new/triangles.png')));
xoff = size(im,2)-fov;
im = im(:,xoff/2+1:end-xoff/2);
bwb = bwboundaries(1-im);

figure(1);clf
imshow(im)
hold on
for i = 1:length(bwb)
    plot(bwb{i}(:,2),bwb{i}(:,1),'r.')
end
line([1 1]*size(im,2)/2,ylim);