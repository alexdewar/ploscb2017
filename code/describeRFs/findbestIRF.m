clear

% % excitatory
% area_ex = 1000; % deg^2
% flat_ex = 0.5;
% cent_ex(1) = 0; % deg
% cent_ex(2) = 0; % deg
% th_ex = pi/2; % rad

initlr = 0.01;

% inhibitory
propr_in = 0;
ecc_in = 0;
cent_in(1) = 0;
cent_in(2) = 0;
th_in = pi/4;

% null
propr_0 = 0;
ecc_0 = 0;
cent_0(1) = 0;
cent_0(2) = 0;
th_0 = 0;

load('vf_kernels.mat','vf_avkernels_r2');
for kern = vf_avkernels_r2
% kern = vf_avkernels_r2(3);

im_size = size(kern.k);

selex = kern.k>0;
rp = regionprops(selex,'Area','Centroid','Eccentricity','Orientation');
area_ex = rp.Area;
cent_ex = rp.Centroid;
ecc_ex = rp.Eccentricity;
th_ex = (pi/180)*rp.Orientation;

irf=makeIRF(im_size, ...
                area_ex,  ecc_ex, cent_ex, th_ex, ...
                propr_in, ecc_in, cent_in, th_in, ...
                propr_0,  ecc_0,  cent_0,  th_0);

figure(2);clf
subplot(1,2,1);
showkernel(kern);
axis square

subplot(1,2,2);
showkernel(irf);
axis square
% hold on
% plot(cent_ex(1),0,'k+')

drawnow
pause(1)

end