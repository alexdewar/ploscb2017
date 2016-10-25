function genellblobs_ablate
    %% constants
    debug = false;

    maxinputs = 100;

    starti = 1;
    nvar = 1000;
    nuvar = 10;
    nwave = 2;
    maxfreq = 30;
    maxamp = 1;
    zrrng = [0.05 0.125];
    av_area = 30;

    suffix = '_ablate';
    a2b = 0.5;

    scale = rand(nvar,1);

    fov = [-135 135; -60 60];
    im_size = [120 270];

    %% ellblob params
    orients = linspace(0,90,nuvar);
    majoraxis = linspace(0,60,nuvar+1);
    majoraxis = majoraxis(2:end);
    el = linspace(fov(2,1),fov(2,2),nuvar);
    %     az = linspace(fov(1,1),fov(1,2),nuvar);

    %% load and tweak kernels
    % load
    load('vf_kernels.mat','vf_avkernels*');
    kerns = [vf_avkernels_r2,vf_avkernels_r4];

    % find mean az (left)
    lk = kerns(cell2mat({kerns.isleft}));
    cent = mean( cell2mat({lk.cent}') ./ ...
                 cell2mat(cellfun(@(x)[size(x,2),size(x,1)],{lk.k},'UniformOutput',false)') );
    az = round(im_size(2)*(cent(1)-0.5));
    
    % resize kernels
    rkerns = resizekernel(kerns,im_size,0.25);
    kcents = cell2mat({kerns.cent}');
    kcents = kcents(:,[2 1]);
    ksz = cell2mat(cellfun(@size,{kerns.k},'UniformOutput',false)');
    rcents = 1+bsxfun(@times,(kcents-1)./(ksz-1),im_size);
    
    r2_ind = 1:length(vf_avkernels_r2);
    r4_ind = length(vf_avkernels_r2)+1:length(kerns);
    
%     figure(1);clf
%     subplot(1,2,1)
%     lefts = cell2mat({vf_avkernels_r4.isleft});
%     showkernels(vf_avkernels_r4(lefts));
%     subplot(1,2,2)
%     showkernels(vf_avkernels_r4(~lefts));
%     keyboard
    
    % gen extra kernels
    [r2ex,r2exc,r2exoi] = genextrakerns(rkerns(:,:,r2_ind),rcents(r2_ind,:),maxinputs,fov);
    [r4ex,r4exc,r4exoi] = genextrakerns(rkerns(:,:,r4_ind),rcents(r4_ind,:),maxinputs,fov);
    [r24ex,r24exc,r24exoi] = genextrakerns(rkerns,rcents,maxinputs,fov);
    
    r2ex_ind = length(kerns)+(1:length(r2exoi));
    r4ex_ind = r2ex_ind(end)+(1:length(r4exoi));
    r24ex_ind = r4ex_ind(end)+(1:length(r24exoi));
    rkerns(:,:,r2ex_ind) = r2ex;
    rkerns(:,:,r4ex_ind) = r4ex;
    rkerns(:,:,r24ex_ind) = r24ex;
    
    %% main loop
    [orients,majoraxis,el,az] = ndgrid(orients,majoraxis,el,az);
    orients = orients(:); majoraxis = majoraxis(:);
    el = el(:); az = az(:);

    minoraxis = majoraxis.*a2b;
    areas = pi.*majoraxis.*minoraxis;
    
    [truearea,aerr] = deal(NaN(nvar,1));
    [amp,freq,phi] = deal(NaN(nvar,nwave));
    trueimsz = prod(im_size./av_area);
    x = NaN([nvar,trueimsz+size(rkerns,3)]);
    if ~debug
        startprogbar(30,nvar);
    end
    i = 1;
    while i <= nvar
        while true
            camp = rand(1,nwave);
            amp(i,:) = range(zrrng)*camp/sum(camp);
            freq(i,:) = randi(maxfreq,1,nwave);
            phi(i,:) = rand(1,nwave)*2*pi;

            blob = bwtrim(~ellblob(amp(i,:),1,freq(i,:),phi(i,:),majoraxis(i),minoraxis(i),pi*orients(i)/180,im_size));

            rp = regionprops(blob,'Orientation','Area','MajorAxisLength','MinorAxisLength');
            aerr(i) = 1-areas(i)./(pi*rp.MajorAxisLength*rp.MinorAxisLength);
    %             fprintf('a1: %f; b1: %f; rat: %f\na2: %f; b2: %f; rat: %f\ndiff: %f\n\n', ...
    %                 majoraxis(i),minoraxis(i),minoraxis(i)/majoraxis(i),rp.MajorAxisLength,rp.MinorAxisLength,rp.MinorAxisLength/rp.MajorAxisLength,aerr(i));

    %             if debug
    %                 figure(1);clf
    %                 imshow(blob);
    %                 keyboard
    %             end

            if length(rp)==1
                break;
            end
        end
    %         if debug
    %             figure(1);clf
    %             imshow(blob);
    %             keyboard
    %         end

        [ys,xs] = find(blob);
        com = [mean(ys) mean(xs)];
        ys = ys-com(1);
        xs = xs-com(2);
        im = false(im_size);
        try
            ind = sub2ind(im_size,min(im_size(1),round(im_size(1)/2+ys)),round(im_size(2)/2+xs));
        catch
            disp('sub2ind')
            continue
        end
        im(ind) = true;

    %         if orient(i)~=0
    %             disp(orient(i))
    %             debug = true;
    %         end
        if debug
            figure(1);clf

            alsubplot(3,1,2,1)
            imshow(im)
        end
        im = imrotate(im,orients(i)-rp.Orientation,'crop');
        im = shiftmat(im,-az(i),round(el(i)));

        truearea(i) = sum(im(:));

        if debug
    %             figure(1);clf

            [xs,ys] = meshgrid(linspace(-im_size(2)/2,im_size(2)/2,im_size(2)),linspace(-im_size(1)/2,im_size(1)/2,im_size(1)));
            sn = sind(orients(i));
            cs = cosd(orients(i));
            cy = ys-el(i);
            ell = hypot((xs*cs-cy*sn)/majoraxis(i),(xs*sn+cy*cs)/minoraxis(i)) <= 1;

            alsubplot(1,1)
            imshow(ell);

            alsubplot(3,1)
            imshow(im);

            keyboard
        end

        im = ~im;

        x(i,trueimsz+1:end) = getacts(im,rkerns);

        rim = blockproc(im,av_area*[1 1],@(x)mean2(x.data));
        x(i,1:trueimsz) = rim(:);

        if ~debug && progbar
            return;
        end
        
        i = i+1;
    end
    % end
    %     im = ~ellblob(scale(i),amp(i,:),freq(i,:),phi(i,:),majoraxis,minoraxis,thoff(i),im_size);
    % 
    
    %% save
    if ~debug
        fname = sprintf('%s/ellblob%s%s.mat',mfiledir,suffix,getenv('JOB_ID'));
        fprintf('Saving to %s...\n',fname);
        save(fname,'x','orients','areas','el','az','majoraxis','minoraxis', ...
             'im_size','scale','amp','phi','freq','maxfreq','maxamp','nvar', ...
             'nwave','starti','maxinputs','r2exc','r2exoi','r4exc','r4exoi', ...
             'r24exc','r24exoi','r2_ind','r4_ind','r2ex_ind','r4ex_ind','r24ex_ind','trueimsz');
    end
end

function [exkerns,excents,exorig]=genextrakerns(rkerns,rcents,maxinputs,fov)
    offrad = 10;

    nex = maxinputs-size(rkerns,3);
    exkerns = NaN([size(rkerns,1),size(rkerns,2),nex]);
    
    rsz = [size(rkerns,1),size(rkerns,2)*360/range(fov(1,:))];
    xd = (rsz(2)-size(rkerns,2))/2;
    zpd = zeros(rsz(1),xd,size(rkerns,3));
    
    canplace = [false(rsz(1),xd), true(rsz(1),size(rkerns,2)), false(rsz(1),xd)];
    
    rkerns = [zpd, rkerns, zpd];
    rcents = bsxfun(@plus,rcents,[0 xd-1]);
    
    [yi,xi] = ndgrid(1:rsz(1),1:rsz(2));
    for i = 1:size(rcents,1)
        exrad(rcents(i,:));
    end
    
    function exrad(c)
        canplace = canplace & hypot(yi-c(1),xi-c(2)) > offrad;
    end

%     figure(1);clf
%     imshow(canplace);
%     keyboard
    
    excents = NaN(nex,2);
    exorig = NaN(1,nex);
    i = 1;
    while i <= nex
        exorig(i) = randi(size(rkerns,3));
        
        crk = rkerns(:,:,exorig(i));
        crc = rcents(exorig(i),:);
        
        placei = find(canplace);
        cplacei = placei(randi(length(placei)));
        [cy,cx] = ind2sub(rsz,cplacei);
        
        ck = circshift(shiftmat(crk,0,cy-crc(1)),[0 round(cx-crc(2))]);
        
%         figure(1);clf
%         imshow(canplace)
%         hold on
%         plot(cx,cy,'g+');
%         pause(.3)
%         drawnow

        % renormalise
        ck = ck(:,xd+1:end-xd);
        ck = sign(ck);
        pos = ck==1;
        neg = ck==-1;
        
        if ~any(pos(:)) || ~any(neg(:))
%             disp(i)
            continue;
        end
        
        ck(pos) = 1/sum(pos(:));
        ck(neg) = -1/sum(neg(:));
        
        exrad([cy cx]);
        
%         figure(1);clf
%         subplot(1,2,1)
%         showkernel(crk,crc([2 1]));
%         axis square
%         subplot(1,2,2)
%         showkernel(ck,[cx-xd cy]);
%         axis square
%         keyboard
        
        exkerns(:,:,i) = ck;
        excents(i,:) = [cy cx-xd];
        
        i = i+1;
    end
end