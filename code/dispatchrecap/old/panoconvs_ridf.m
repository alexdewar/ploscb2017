function panoconvs_ridf(dosave)
    if ~nargin
        dosave = false;
    end
    
    doall = false;
    ymax = 0.6*ones(1,4);
    xtick = [-180 -90 0 90 180];
    if doall
        sz = [7 10];
    else
        sz = [7 5];
    end

    dname = 'antoinestim';
    
    load('vf_kernels','vf_avkernels*');

    fulldname = fullfile(mfiledir,dname); %,'touse'
    if ~doall
        fulldname = [fulldname,'/touse'];
    end
    d = [dir(fullfile(fulldname,'*.jpg'));dir(fullfile(fulldname,'*.png'))];

    for i = 1:length(d)
        fname = fullfile(fulldname,d(i).name);
        fnnoext = d(i).name(1:end-4);
        
        im = imread(fname);
        if size (im,3)>1;
            im = rgb2gray(im);
        end
        im = im2double(im);

        figure(i);clf
        
        if doall
            alsubplot(2,1,1,1)
            imagesc(im);
            colormap gray
            title(fnnoext,'Interpreter','none')
            axis off

            alsubplot(2,1)
        end
        
        [rr2,thr2]=vf_ridf(im,vf_avkernels_r2);
        [rr4,thr4]=vf_ridf(im,vf_avkernels_r4);
        
        set(gca,'FontSize',8);
        hold on
                
        plot(thr2,rr2,thr4,rr4);
        
        if doall
            cymax = max([rr2,rr4]);
        else
            cymax = ymax(i);
        end
        line([0 0],[0 cymax],'Color','k','LineStyle','--');
        line([-90 -90],[0 cymax],'Color','k','LineStyle','--');
        
        diffr2 = rr2(thr2==-90);
        diffr4 = rr4(thr4==-90);
        
        set(gca,'XTick',xtick);
        
        axis tight
        ylabel('r.m.s. difference')
%         xlabel('Orientation (\circ)')

        format_ticks(gca);

        if dosave
            savefig(sprintf('%s_%s_r2_%.5f_r4_%.5f',dname,fnnoext,diffr2,diffr4),sz)
        end
    end
    
    if dosave
        close all
    end
end

function [ridf,ths]=vf_ridf(im,kerns)
    kerns = cell2mat(shiftdim({kerns.k},-1));
    
    ksz = size(kerns);
    rim = imresize(im,[ksz(1),ceil(360*ksz(2)/270)]);
    xoff = (size(rim,2)-ksz(2))/2;
    acts = NaN(ksz(3),size(rim,2));
    for i = 1:size(rim,2)
        crim = circshift(rim,[0 i-1-ceil(size(rim,2)/2)]);
        
        acts(:,i) = shiftdim(sum(sum(bsxfun(@times,crim(:,xoff+(1:ksz(2))),kerns)),2));
    end
    
%     for i = 1:size(im,2)
%         [acts(:,i),kerns] = getneuronactivations(circshift(im,[0 i-1-ceil(size(im,2)/2)]),kerns);
%     end
    
    ridf = sqrt(mean(bsxfun(@minus,acts,acts(:,1+size(rim,2)/2)).^2));
    ridf = [ridf,ridf(1)];
    
    ths = linspace(-180,180,size(rim,2)+1);
end