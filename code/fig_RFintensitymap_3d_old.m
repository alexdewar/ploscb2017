function fig_RFintensitymap_3d
    r2_or_r4 = 'r4';

    load('vf_kernels.mat');
    main(eval(['vf_avkernels_' r2_or_r4]));
end

function main(kstruct)
    xlims = [-135 135];
    ylims = [-60 60];
    filtsize = 3;
    imscale = 1;
    
    cylht = 2;
    cylr  = 0.5;

    %% make intensity map image
    kstruct = kstruct(cell2mat({kstruct.isleft})); % & cell2mat({kstruct.glomnum})==5);
    kernels = {kstruct.k};
    imsz = size(kstruct(1).k);
    
    mapim = ones([imscale.*imsz,3]);
    ecol = [1 0 0];
    icol = [0 0 1];
    nkern = numel(kernels);
    for i = 1:nkern
        curk = kernels{i};
        addtoim(curk>0,ecol);
        addtoim(curk<0,icol);
    end
    
    function addtoim(cim,col)
%         figure(5);imshow(cim)
%         cim = medfilt2(imresize(cim,imscale),filtsize*[1 1]);

        bnd = bwboundaries(medfilt2(cim,filtsize*[1 1]));
        outlines = false(size(cim));
        for b = 1:numel(bnd)
            outlines(sub2ind(size(cim),bnd{b}(:,1),bnd{b}(:,2))) = true;
        end
        
        fim = imfill(outlines,'holes');
        newim = ones([size(cim),3]);
        alpha = zeros(size(cim));
        for r = 1:size(newim,1)
            newim(r,fim(r,:),1) = col(1);
            newim(r,fim(r,:),2) = col(2);
            newim(r,fim(r,:),3) = col(3);
            alpha(r,fim(r,:)) = 1/nkern;
            alpha(r,outlines(r,:)) = 1;
        end

        mapim = bsxfun(@times,newim,alpha) + bsxfun(@times,mapim,1-alpha);
    end

    %% convert to polygons
%     mapim(:,end,:) = 0;
    
    vsz = prod(imsz);
    
    cols = reshape(mapim,vsz,3);
    notwhites = ~all(cols==1,2);
    cols = cols(notwhites,:);
    
    [yi,xi] = ind2sub(imsz,1:vsz);
    Xim = bsxfun(@plus,xi(notwhites)',[-1 0 0 -1]);
    Yim = bsxfun(@plus,yi(notwhites)',[-1 -1 0 0]);
    
    %% wrap round cylinder
    xlimsr = pi*xlims/180;
    xrng = range(xlimsr);
    yrng = range(pi*ylims/180);
    radperpix = 1./((imsz-1).*[yrng xrng]);
    pxsz = 2*cylr*tan(radperpix/2);
    
    Z3 = Yim*yrng*cylht./(imsz(1)*2*atan(cylht/(2*cylr)));
    [X3,Y3] = pol2cart(Xim*xrng./imsz(2),cylr);
    
    figure(1);clf
    hold on
    for i = 1:size(Xim,1)
        fill3(X3,Y3,Z3,cols(i,:),'linestyle','none');
    end

    ths = linspace(0,2*pi,100);
    [linex,liney] = pol2cart(ths,cylr);
    linez = zeros(size(linex));

    line(linex,liney,linez);
    line(linex,liney,linez+cylht);
    line([0 0],-cylr*[1 1],[0 cylht]);
    line([0 0],cylr*[1 1],[0 cylht]);
end