function findcompattern(gennew,dosave)
    if nargin < 1
        gennew = false;
    end
    if nargin < 2
        dosave = false;
    end
    
    simblobs = false;
    
    blobalpha = 0.5;

    fov = [120 270];
    
    comdiff = 20;

    load('vf_kernels.mat','vf_avkernels_r2');
    kerns = vf_avkernels_r2;
    ks = cell2mat(shiftdim({kerns.k},-1));
    imsz = [size(ks,1),size(ks,2)];
    bigimsz = 10*imsz;
%     blobimsz = [imsz(1), imsz(2)*(360/fov(2))];
%     xoff = (blobimsz(2)-imsz(2))/2;
    
    comdiff = comdiff*(imsz(1)/fov(1));
    
    cents = cell2mat({kerns(cell2mat({kerns.isleft})).cent}');
    az = round(mean(cents(1,:)) - imsz(2)/2);
    
    a2b = 0.5;
    nwave = 2;
    
    % scale,majoraxis,thoff,amp,freq,phi
    maxes = [  1 30*(imsz(1)/fov(1)) pi, 0.125 30 2*pi]';
    mins  = [0.1  5*(imsz(1)/fov(1))  0,  0.05  1    0]';
    
    maxes = [maxes(1:3);repmat(maxes(4:end),nwave,1)];
    mins = [mins(1:3);repmat(mins(4:end),nwave,1)];
    rng = maxes-mins;
%     param1 = mins(1:3)+rng(1:3).*rand(3,1);
%     wparam1 = bsxfun(@plus,mins(4:6)',bsxfun(@times,rng(4:6)',rand(nwave,3)));
%     wparam1(:,2) = round(wparam1(:,2));
%     wparam1 = wparam1(:);
    
%     param1 = rand(9,1);
%     wparam1 = bsxfun(@plus,mins(4:6)',bsxfun(@times,rng(4:6)',rand(nwave,3)));
%     wparam1(:,2) = round(wparam1(:,2));
%     wparam1 = wparam1(:);
    
%     param2 = mins(1:3)+rng(1:3).*rand(3,1);
%     wparam2 = bsxfun(@plus,mins(4:6)',bsxfun(@times,rng(4:6)',rand(nwave,3)));
%     wparam2(:,2) = round(wparam2(:,2));
%     wparam2 = wparam2(:);

%     errval=errfunc(x0);
%     keyboard

    datafn = sprintf('%s/fcp.mat',mfiledir);
    if ~gennew && exist(datafn,'file')
        load(datafn);
    else
        disp('searching for new blob combo')
        
        opts = optimset('fminsearch');
        [opts.MaxIter,opts.MaxFunEvals] = deal('2000*numberofvariables');
        x0 = rand(6+6*nwave,1);
        [xout,fval] = fminsearch(@(p)errfunc(p,imsz,simblobs),x0)

        [~,evact,evparam] = errfunc(xout);

        [~,bacts1,blob1] = blobacts(xout(1:length(xout)/2),0,bigimsz);
        [~,bacts2,blob2] = blobacts(xout(length(xout)/2+1:end),comdiff,bigimsz);
        
        save(datafn,'xout','fval','evact','evparam','blob1','blob2','bacts1','bacts2');
    end
    
    y = linspace(-fov(1)/2,fov(1)/2,imsz(1));
    x = linspace(-fov(2)/2,fov(2)/2,imsz(2));

    blob1 = im2double(blob1);
    blob1(blob1==1) = NaN;
    blob2 = im2double(blob2);
    blob2(blob2==1) = NaN;
    blob2 = repmat(blob2,[1 1 3]);
    blob2(:,:,2) = 1-blob2(:,:,2);
%     blob2 = repmat(blob2,[1 1 3]);
%     blob2(:,:,4) = blobalpha*~blob2(:,:,1);
%     blob2(:,:,[1 3]) = false;
    
    blobs = imagealpha(blob1,blobalpha,blob2,blobalpha);
    
    figure(1);clf
    subplot(2,1,1)
    image(x,y,blobs);
    colormap gray
    axis equal
    
    subplot(2,1,2)
    h=bar([bacts1, bacts2]);
    c1 = 1-blobalpha*(1-[0 0 0]);
    c2 = 1-blobalpha*(1-[0 1 0]);
    set(h(1),'FaceColor',c1);
    set(h(2),'FaceColor',c2);
    ylim([-1 1]);
    
    xlabs = cell(size(kerns));
    for i = 1:length(kerns)
        if kerns(i).isleft
            xlabs{i} = sprintf('L%d',kerns(i).glomnum);
        else
            xlabs{i} = sprintf('R%d',kerns(i).glomnum);
        end
    end
    set(gca,'XTick',1:length(kerns),'XTickLabel',xlabs);
    rdiff = getRMSdiff(bacts1,bacts2);
    text(length(kerns)+0.5,1,sprintf('mean difference = %.1f%%',100*rdiff), ...
         'HorizontalAlignment','right','VerticalAlignment','top');
    
%     image(x,y,255*blob1);
%     colormap gray
%     axis equal
%     
%     subplot(1,3,3)
%     image(x,y,255*blob2)
%     colormap gray
%     axis equal
    
%     subplot(1,2,2)
%     image(x,y,255*blob2);
%     axis equal
%     colormap gray
    
%     subplot(2,-1,2)
%     image(x,y,255*blob2);
%     colormap gray
    
%     fprintf('difference: %f\n\n',getRMSdiff(bacts1,bacts2));
    
%     subplot(4,1,3)
%     imshow(abs(blob1-blob2))
%     
% %     subplot(4,1,4)
%     showkernels(kerns,[],.1,fov(2)*[-0.5 0.5],fov(1)*[-0.5 0.5]);
    
%     keyboard

    dump2base(true)
    
    if dosave
        savefig('fcp',[20 10]);
    end
    
    function [errval,evact,evparam]=errfunc(params,cimsz,simordissim)
        [~,acts1] = blobacts(params(1:length(params)/2),0,cimsz);
        [~,acts2] = blobacts(params(length(params)/2+1:end),~simordissim*comdiff,cimsz);
        evact = mean(abs((acts1-acts2)/2));

%         evim  = mean(abs(im1(:)-im2(:)));
        evparam = mean(abs(params(1:length(params)/2)-params(length(params)/2+1:end)));
%         disp(evact)
        if simordissim
            errval = evparam/evact;
        else
            errval = evact/evparam;
        end
%         disp(errval);
%         disp(params)
    end

    function [im,acts,trueim]=blobacts(b_x,vcomoff,cimsz)
        [im,trueim] = makeblob(b_x,vcomoff,cimsz);
        acts = getacts(trueim,ks);
    end
    
    function [im,trueim]=makeblob(b_x,vcomoff,cimsz)
        b_x = max(mins,min(maxes,rng.*(b_x+mins)));
        
        b_param = b_x(1:3);
        b_wparam = b_x(3+(1:nwave*3));
        b_wparam = reshape(b_wparam,length(b_wparam)/nwave,nwave);
%         for w = 1:size(b_wparam,2)
%             b_wparam(:,w) = max(mins(3+w),min(maxes(3+w), ...
%                                 b_x(3+(w-1)*nwave+(1:nwave))));
% %             b_wparam(:,w) = b_x(3+(w-1)*nwave+(1:nwave));
%         end
        
        % scale,majoraxis,thoff,amp,freq,phi
        % (scale,amp,freq,phi,majoraxis,minoraxis,thoff,im_size)
        im = ellblob(b_param(1),b_wparam(1,:),b_wparam(2,:),b_wparam(3,:), ...
                     b_param(2),b_param(2)*a2b,b_param(3),cimsz);
        if all(im(:))
            error('no blob!')
        end
        
        [ys,~] = find(~im);
        cvcom = mean(ys);
        
        im = circshift(im,round(cvcom-cimsz(1)/2-vcomoff));
        
        bwid = cimsz(2)*fov(2)/(2*360);
        trueim = circshift(im,[0,az]);
        trueim = trueim(:,round((cimsz(2)-bwid)/2+(1:round(bwid))));
        
%         figure(1);clf
%         imshow(im)
%         keyboard
    end

end