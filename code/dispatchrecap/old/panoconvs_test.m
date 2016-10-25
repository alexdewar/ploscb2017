clear

dname = [mfiledir '/antoinestim/touse/'];
fname = '09_2_*.png';
d = dir([dname fname]);
disp(length(d))
im = im2double(rgb2gray(imread([dname d(1).name])));

fov = 270;

imsz2 = [size(im,1),size(im,2)*(fov/360)];

load('vf_kernels.mat','vf_avkernels_r2');
kerns = vf_avkernels_r2;
rkerns = resizekernel(kerns,imsz2,.25);

xoff = (size(im,2)-imsz2(2))/2;
rots = -180:179;
acts = NaN(length(kerns),360);
for i = 1:360
    crim = circshift(im,[0 rots(i)]);
    acts(:,i) = getacts(crim(:,xoff+(1:imsz2(2))),rkerns);
end

%%
figure(1);clf
alsubplot(1,2,1,1)
plot(acts')
alsubplot(1,2)
plot(mean(acts))