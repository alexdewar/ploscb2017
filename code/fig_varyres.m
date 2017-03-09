function fig_varyres
    load('ofstad_views_full');
    load('vf_kernels');
    
    main(views,lr_views,superlr_views,d,xg,yg,vf_avkernels_r2,vf_avkernels_r4,thresh);
end

function main(views,lr_views,superlr_views,d,xg,yg,vf_avkernels_r2,vf_avkernels_r4,thresh)
    fov = 270; % degrees (horizontally)
    catchment_area_thresh = pi/2;
    
    % trim edges of views
    eoff = round(size(lr_views,2)*(360-fov)/(2*360));
%     views = lr_views(:,eoff+1:end-eoff,:);
    
    % get ring neurons' "views"
    kernels = {vf_avkernels_r2(:).k,vf_avkernels_r4(:).k};
    r2_ind = 1:numel(vf_avkernels_r2);
    r4_ind = numel(vf_avkernels_r2)+1:numel(kernels);
    
    rkernels = NaN(size(lr_views,1),size(lr_views,2),numel(kernels));
    for i = 1:numel(kernels)
        rkernels(:,:,i) = resizekernel(kernels{i},[size(lr_views,1),size(lr_views,2)],thresh);
    end
    
    k_views = NaN(numel(kernels),1,size(views,3));
    for i = 1:size(views,3)
        for j = 1:numel(kernels)
            k_views(j,1,i) = sum(sum(lr_views(:,:,i).*rkernels(:,:,j)));
        end
    end
    
    % main loop
    [refx,refy] = pol2cart((pi/4)+[0 pi/2 pi 1.5*pi],d/4);
%     [refx,refy] = meshgrid(linspace(-d/2,d/2,7));
%     selref = hypot(refx(:),refy(:)) <= d;
%     refx = refx(selref);
%     refy = refy(selref);
    mini = NaN(size(refx));
    allidf = NaN(34);
    figure(2);clf
    hold on
    for i = 1:numel(refx)
        % reference view
        [~,xi] = min(abs(refx(i)-xg));
        [~,yi] = min(abs(refy(i)-yg));
        ref_view_i = xg==xg(xi) & yg==yg(yi);

        %% RMS IDF
%         figure(i*10+1);clf
%         idfplot(views,'RMS image difference');

        %% R2 IDF
        [mx,my,cidf,bndi]=idfplot(k_views(r2_ind,:,:),'R2 neurons');
        
        [~,mini(i)] = min(cidf(:));
        
        switch i
            case 1
                sel = mx > 0 & my > 0;
            case 2
                sel = mx <= 0 & my > 0;
            case 3
                sel = mx <= 0 & my <= 0;
            case 4
                sel = mx > 0 & my <= 0;
        end
        allidf(sel) = cidf(sel);
        
%         if i==1 || i==4
%             sel = xg > 0;
%             xg = xg(sel);
%             yg = yg(sel);
%         else
%             sel = xg <= 0;
%             xg = xg(sel);
%             yg = yg(sel);
%         end
%         if i < 3
%             sel = yg > 0;
%             xg = xg(sel);
%             yg = yg(sel);
%         else
%             sel = yg <= 0;
%             xg = xg(sel);
%             yg = yg(sel);
%         end

        %% R4 IDF
%         figure(i*10+3);clf
%         idfplot(k_views(r4_ind,:,:),'R4 neurons');
% 
%         %% R2 & R4 IDF
%         figure(i*10+4);clf
%         idfplot(k_views,'R2 & R4 neurons');
% 
%         %% low-res RMS IDF
%         figure(i*10+5);clf
%         idfplot(superlr_views,'low-res RMS');
    end
    
    contourf(mx,my,allidf);
    colormap gray
    
    plot(mx(mini),my(mini),'r+');
    
    r = d/2;
    line([0 0],[-r r],'LineStyle','--')
    line([-r r],[0 0],'LineStyle','--')
    ths = linspace(0,2*pi,1000);
    [rimx,rimy] = pol2cart(ths,r);
    line(rimx,rimy,'LineStyle','--')
    
    ylim([-r r]);
    xlim([-r r]);
    keyboard
    
    %% plotting function
    function [mx,my,idf,bndi]=idfplot(cviews,ttl)
        rms_diffs = shiftdim(getRMSdiff(cviews,cviews(:,:,ref_view_i)));
        
        [idf,mx,my] = makeim(xg,yg,rms_diffs);

        heads = getIDFheads(xg,yg,rms_diffs);
        [heads_im,mxh,myh] = makeim(xg,yg,heads);
        success = abs(heads_im) < catchment_area_thresh;
        success(ref_view_i) = true;
        success = bwlabeln(success);
        success = success==success(ref_view_i);
        bnd = bwboundaries(success);
        bnd = bnd{:};
        bndi = sub2ind(size(heads_im),bnd(:,1),bnd(:,2));
%         rp = regionprops(success,'Area');
% 
%         hold on
%         surf(mx,my,idf);
%         plot(mxh(bndi),myh(bndi),'r');
%         title(sprintf('%s (Area: %.2f)',ttl,rp.Area));
    end
end