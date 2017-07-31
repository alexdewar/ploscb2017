function barpanoconv_pic(dosave)
if ~nargin
    dosave=false;
end

imfns = {'tree1.jpg','tree2.jpg','tree3.jpg'};

thresh = 0.25;
fov = 270;

dosavefigdat = true;

datadir = fullfile(mfiledir,'../drosodata/barextra');
if ~exist(datadir,'dir')
    mkdir(datadir);
end

for i = 1:length(imfns)

    imfn = imfns{i};
    disp(imfn)

    im = im2double(rgb2gray(imread(imfn)));
    imsz = size(im);

    ksz = [imsz(1), imsz(2)*fov/360];

    load('vf_kernels.mat','vf_avkernels_*');
    rkerns_r2 = resizekernel(vf_avkernels_r2,ksz,thresh);
    r2_l = cell2mat({vf_avkernels_r2.isleft});
    rkerns_r4 = resizekernel(vf_avkernels_r4,ksz,thresh);
    r4_l = cell2mat({vf_avkernels_r4.isleft});

    datafn = sprintf('%s/vw_barpic_%s_%s.mat',datadir,imfn(1:end-4),mfilehash);
    if exist(datafn,'file')
        load(datafn);
    else
        pfunc = @panoconv_all;

        disp('R2 L')
        vals_r2_l = pfunc(im,rkerns_r2(:,:,r2_l),fov);
        disp('R2 R')
        vals_r2_r = pfunc(im,rkerns_r2(:,:,~r2_l),fov);
        disp('R4 L')
        vals_r4_l = pfunc(im,rkerns_r4(:,:,r4_l),fov);
        disp('R4 R')
        [vals_r4_r,ths] = pfunc(im,rkerns_r4(:,:,~r4_l),fov);
        disp('.')

        if dosavefigdat
            save(datafn,'vals*','ths');
        end
    end

    %% endhash
    
    ths(end+1) = -ths(1);
    vals_r2_l(:,end+1) = vals_r2_l(:,1);
    vals_r2_r(:,end+1) = vals_r2_r(:,1);
    vals_r4_l(:,end+1) = vals_r4_l(:,1);
    vals_r4_r(:,end+1) = vals_r4_r(:,1);

    figure(i);clf
    
    subplot(5,1,1)
    imagesc(im)
    axis off
    
    subplot(5,1,2)
    indivplot(ths,vals_r2_l,'R2 L');
    
    subplot(5,1,3)
    indivplot(ths,vals_r2_r,'R2 R');
    
    subplot(5,1,4)
    indivplot(ths,vals_r4_l,'R4 L');
    
    subplot(5,1,5)
    indivplot(ths,vals_r4_r,'R4 R');

    if dosave
        savefig(sprintf('vw_barpic_%s',imfn(1:end-4)),[20 6]);
    end

end

function indivplot(ths,vals,label)
hold on
plot(ths,vals)
plot(ths,mean(vals),'k--','LineWidth',3);
xlim([ths(1) ths(end)])
ylim([-0.5 0.5])
set(gca,'FontSize',8,'XTick',-180:45:180);
title(label)