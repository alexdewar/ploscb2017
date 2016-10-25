function panoconvs_ridf_scores(dosave)
    if ~nargin
        dosave = false;
    end
    
    doall = false;
    showprogbar = isempty(dbstatus);
    
    ncol = 2;
    ntick = 3;
    tickmax = 0.4;
    
%     ymax = 0.6*ones(1,4);
%     xtick = [-180 -90 0 90 180];
%     if doall
        sz = [16 7];
%     else
%         sz = [7 5];
%     end
    
    xlo = 0.05;
    showrng = [45 225];
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

    [diffr2,stdr2,diffr4,stdr4] = deal(NaN(ceil(length(d)/ncol),ncol));
    if showprogbar
        startprogbar(1,length(d));
    end
    
    crng = round(showrng*im_size(2)/360);
    
    figure(1);clf
    for i = 1:ncol
        subplot(1,ncol,i);
        set(gca,'FontSize',8);
        hold on
        
        ims = ones(im_size(1)*size(diffr2,1),range(crng)+1);
        for j = 1:size(diffr2,1)
            ind = (i-1)*size(diffr2,1)+j;
            if ind>numel(fname)
                break;
            end
            
            im = imread(fullfile(fulldname,fname{ind}));
            if size (im,3)>1;
                im = rgb2gray(im);
            end
            im = im2double(im);

            ims((j-1)*im_size(1)+(1:im_size(1)),:) = im(:,1+(crng(1):crng(2)));

            [rr2,thr2]=vf_ridf(im,vf_avkernels_r2);
            diffr2(j,i) = rr2(thr2==-90);
            stdr2(j,i) = std(rr2)/sqrt(length(rr2));
            
            [rr4,thr4]=vf_ridf(im,vf_avkernels_r4);
            diffr4(j,i) = rr4(thr4==-90);
            stdr4(j,i) = std(rr4)/sqrt(length(rr4));

            if showprogbar && progbar
                return;
            end
        end
        
        horzerr(-diffr2(:,i),-stdr2(:,i),-xlo)
        horzerr(diffr4(:,i),stdr4(:,i),xlo)
        
        x = linspace(-xlo,xlo,im_size(2));
        y = 0.5+linspace(size(diffr2,1),0,size(ims,1));
        image(x,y,ims*255);
        colormap gray

%         set(gca,'YTick',[]);
        ticks = [linspace(-tickmax,-xlo,ntick),linspace(xlo,tickmax,ntick)];
        set(gca,'XTick',ticks,'XTickLabel',num2str(abs(ticks')-xlo),'YTick',[]);
%         xtick = get(gca,'XTick');
%         set(gca,'XTick',xtick(xtick>=0));
        xlim([-1 1] * (xlo+1.05*max([diffr2(:)+stdr2(:);diffr4(:)+stdr4(:)])));
        xlabel('r.m.s. difference between 0^\circ and 90^\circ')
    end
    
    if dosave
        savefig('pattern_ridf',sz);
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

function horzerr(x,h,xlo)
    x = x+xlo;
    y = size(x,1):-1:1;
    barh(y,x,'w');
    for i = 1:numel(x)
        line(x(i)+[0 h(i)],[y(i) y(i)],'Color','k')
        line(x(i)+h(i)*[1 1],y(i)+[-0.2 0.2],'Color','k')
    end
end