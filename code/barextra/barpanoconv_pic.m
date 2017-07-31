function barpanoconv_pic(dosave)
if ~nargin
    dosave=false;
end

imfns = {'tree1.jpg','tree2.jpg'};

thresh = 0.25;
fov = 270;

dosavefigdat = true;

datadir = fullfile(mfiledir,'../drosodata/barextra');
if ~exist(datadir,'dir')
    mkdir(datadir);
end

[ims,datafn] = deal(cell(1,length(imfns)));

for i = 1:length(imfns)
    imfn = imfns{i};
    disp(imfn)
    
    im = im2double(rgb2gray(imread(imfn)));
    ims{i} = im;
    
    datafn{i} = sprintf('%s/vw_barpic_%s_%s.mat',datadir,imfn(1:end-4),mfilehash);
    if ~exist(datafn{i},'file')
        
        imsz = size(im);
        ksz = [imsz(1), imsz(2)*fov/360];
        
        load('vf_kernels.mat','vf_avkernels_r4');
        rkerns_r4 = resizekernel(vf_avkernels_r4,ksz,thresh);
        
        disp('R4')
        [vals_r4,ths] = panoconv_all(im,rkerns_r4,fov);
        disp('.')
        
        if dosavefigdat
            save(datafn{i},'vals_*','ths');
        end
    end
    
    %% endhash
end

mvals = cell(1,length(imfns));
maxact = 0;
for i = 1:length(imfns)
    load(datafn{i});
    
    ths(end+1) = -ths(1);
    vals_r4(:,end+1) = vals_r4(:,1);
    
    nkern = size(vals_r4,1);
    cmvals = [mean(vals_r4(1:nkern/2,:)); mean(vals_r4(nkern/2+1:end,:))];
    maxact = max(maxact,max(abs(cmvals(:))));
    
    mvals{i} = cmvals;
end

figure(1);clf
for i = 1:length(imfns)
    cvals = mvals{i} / maxact;
    
    subplot(length(imfns),1,i)
    hold on
    
    imagesc(ths([1 end]),[-1 1],flipud(ims{i}));
    colormap gray
    axis tight
    pbaspect([size(ims{i},2),size(ims{i},1),1])
    
    plot(ths,[mean(cvals); diff(cvals)],'LineWidth',3);
    plot(ths,cvals,'LineWidth',2,'LineStyle','--');

    set(gca,'XTick',-180:45:180);
    xlabel('Angle')
    ylabel('Mean activation')
    
    legend({'mean', 'diff', 'left', 'right'})
    
    title(imfns{i})
end

if dosave
    savefig('vw_barpic',[20 15]);
end