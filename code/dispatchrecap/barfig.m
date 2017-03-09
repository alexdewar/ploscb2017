function barfig(dosave)
    if ~nargin
        dosave = false;
    end
    
    tgain = -1;
    
    ksz = [120 270];
    
    load('vf_kernels','vf_avkernels*');
    r2k = resizekernel(vf_avkernels_r2(cell2mat({vf_avkernels_r2.isleft})),ksz,.25);
    r4k = resizekernel(vf_avkernels_r4(cell2mat({vf_avkernels_r4.isleft})),ksz,.25);
    
    barfs = dir(fullfile(mfiledir,'../barextra/vw_bar*.png'));

    for i = 1:length(barfs)
        im = im2double(imread(fullfile(mfiledir,'../barextra',barfs(i).name)));
        
        figure(i);clf
        set(gca,'FontSize',8);
        alsubplot(4,1,1,1);
        imagesc(linspace(-180,180,size(im,2)),linspace(-60,60,size(im,1)),im);
        colormap gray
        
        actsr2 = getracts(im,r2k);
        actsr4 = getracts(im,r4k);
        
        alsubplot(2,1)
        plotacts(mean(actsr2));
        
        alsubplot(3,1)
        plotacts(mean(actsr4));
        
        hiactsr2 = mean(max(actsr2,0));
        hiactsr4 = mean(max(actsr4,0));
        alsubplot(4,1)
        plot(linspace(-180,180,length(actsr2)),tgain*180*(hiactsr2-hiactsr2(end:-1:1)+hiactsr4-hiactsr4(end:-1:1)));
        xlim([-180 180])
        ylim([-180 180])
        set(gca,'XTick',-180:90:180,'YTick',-180:90:180)
        
        if dosave
            savefig(barfs(i).name(1:end-4),[7 5]);
        end
    end
end

function plotacts(macts)
    yl = [-.25 .25];
    
    ths = linspace(-180,180,length(macts));
    
    plot(ths,macts,'k',ths,macts(end:-1:1),'k--');
%     ylim(yl);
    xlim([-180 180]);
    box off
    
    set(gca,'XTick',-180:90:180);
%     format_ticks(gca);
end

function racts=getracts(im,rkerns)
    racts = NaN(size(rkerns,3),size(im,2));
    for i = 1:size(im,2)
        racts(:,i) = getacts(cshiftcut(im,size(rkerns,2),1-i),rkerns);
    end
    
%     macts = mean(acts);
%     macts = [macts,macts(1)];
end