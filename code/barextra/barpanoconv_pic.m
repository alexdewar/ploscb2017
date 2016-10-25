function barpanoconv_pic(dosave)
if ~nargin
    dosave=false;
end

imfn = 'panotrees.jpg';

thresh = 0.25;
fov = 270;

im = im2double(rgb2gray(imread(imfn)));
imsz = size(im);

ksz = [imsz(1), imsz(2)*fov/360];

load('vf_kernels.mat','vf_avkernels_*');
rkerns_r2 = resizekernel(vf_avkernels_r2,ksz,thresh);
r2_l = cell2mat({vf_avkernels_r2.isleft});
rkerns_r4 = resizekernel(vf_avkernels_r4,ksz,thresh);
r4_l = cell2mat({vf_avkernels_r4.isleft});

datafn = sprintf('%s/vw_barpic_%s.mat',mfiledir,mfilehash);
if exist(datafn,'file')
    load(datafn);
else
    pfunc = @panoconv;
    
    disp('R2 L')
    vals_r2_l = pfunc(im,rkerns_r2,fov);
    disp('R2 R')
    vals_r2_r = pfunc(im,rkerns_r2(:,:,~r2_l),fov);
    disp('R4 L')
    vals_r4_l = pfunc(im,rkerns_r4,fov);
    disp('R4 R')
    [vals_r4_r,ths] = pfunc(im,rkerns_r4(:,:,~r4_l),fov);
    disp('.')

    save(datafn,'vals*','ths');
end

%% endhash

figure(1);clf
% plot(ths,vals_r2_l,ths,vals_r2_r,'b--', ...
%      ths,vals_r4_l,'g',ths,vals_r4_r,'g--')
hold on
plot(ths,vals_r2_l,ths,vals_r2_r,'b--')
plot(ths,vals_r4_l,'g',ths,vals_r4_r,'g--')
% alpha 0.5

xlim([-180 180])

% axis square

% ylim(0.5*[-1 1])

set(gca,'FontSize',8,'XTick',-180:180:180); %,'YTick',-0.5:0.5:0.5);
% plot(ths,mean(vals(lefts,:)),'k',ths,mean(vals(~lefts,:)),'k--')


% subplot(1,2,2)
% plot(ths,mean(vals(~lefts,:)),'k--');

% ylabel('Activation')
% axis tight

if dosave
    savefig('vw_barpic',[20 10])
end