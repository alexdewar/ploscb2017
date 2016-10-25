function panoconvs_ridf_scores_spp(dosave)
    if ~nargin
        dosave = false;
    end
    
    doall = false;
    showprogbar = false;
    
    nrow = 3;
    ntick = 3;
    tickmax = 0.4;
    
    ytick = 0:0.2:1;
    
%     ymax = 0.6*ones(1,4);
%     xtick = [-180 -90 0 90 180];
%     if doall
        sz = [16 30];
%     else
%         sz = [7 5];
%     end
    
    ylo = 0.1;
    showrng = [135 315];
    im_size = [120 360];

    dname = 'antoinestim';
    
    load('vf_kernels','vf_avkernels*');

    fulldname = fullfile(mfiledir,dname);
    if ~doall
        fulldname = [fulldname,'/new'];
    end
    d = [dir(fullfile(fulldname,'*.jpg'));dir(fullfile(fulldname,'*.png'))];
    fname = sort({d.name});
	fname = fname(end:-1:1);

    spp = NaN(ceil(length(d)/nrow),2,nrow);
    if showprogbar
        startprogbar(1,length(d));
    end
    
    crng = round(showrng*im_size(2)/360);
    
    csz = range(crng)+1;
    figure(1);clf
    for i = 1:nrow
        subplot(nrow,1,i);
        set(gca,'FontSize',8);
        hold on
        ims = ones(im_size(1),csz*size(spp,1));
        for j = 1:size(spp,1)
            ind = (i-1)*size(spp,1)+j;
            if ind>numel(fname)
                break;
            end
            
            im = imread(fullfile(fulldname,fname{ind}));
            if size (im,3)>1;
                im = rgb2gray(im);
            end
            im = im2double(im);

            ims(:,(j-1)*csz+(1:csz)) = im(:,1+(crng(1):crng(2)));

%             [rr2,thr2]=vf_ridf(im,vf_avkernels_r2);
%             diffr2(j,i) = rr2(thr2==-90);
%             stdr2(j,i) = std(rr2)/sqrt(length(rr2));
            
            spp(j,1,i)=vf_spp(im,vf_avkernels_r4);
            spp(j,2,i)=vf_spp(im,vf_avkernels_r2);

            if showprogbar && progbar
                return;
            end
        end
        
        bar(ylo*sign(spp(:,:,i))+spp(:,:,i),'w');
        
        x = 0.5+linspace(0,size(spp,1),size(ims,2));
        y = linspace(ylo,-ylo,size(ims,1));
        image(x,y,ims*255);
        colormap gray
        
        yt = [-ytick(end:-1:1)-ylo,ytick+ylo];
        ytstr = num2str([ytick(end:-1:1),ytick]');
        set(gca,'XTick',[],'YTick',yt,'YTickLabel',ytstr);
        
    end
    
    if dosave
        savefig('pattern_spp',sz);
    end
end

function spp=vf_spp(im,kerns)
    kerns = cell2mat(shiftdim({kerns.k},-1));
    
    ksz = size(kerns);
    rim = imresize(im,[ksz(1),ceil(360*ksz(2)/270)]);
    xoff = (size(rim,2)-ksz(2))/2;
    ind = xoff+1:size(rim,2)-xoff;
    acts0 = mean(getacts(rim(:,ind),kerns));
    crim = circshift(rim,[0 round(size(rim,2)/4)]);
    acts90 = mean(getacts(crim(:,ind),kerns));
    
%     disp([acts0,acts90]);
    spp = (acts0-acts90)./(abs(acts0)+abs(acts90));
end
