function im=getviewfast(x,y,z,th,X,Y,Z,av_area,el_max,orig_im_size,pitch)
% function im=getviewfast(x,y,z,th,X,Y,Z,av_area,el_max,orig_im_size,pitch)
%     fillthresh = 0.2;

    if nargin < 10 || isempty(orig_im_size)
        orig_im_size = [170 900];
    end
    if nargin < 9 || isempty(el_max)
        el_max = 5*pi/12;
    else
        el_max = pi*el_max/180;
    end

    if nargin < 11
        pitch = el_max/2;
    else
        pitch = pi*pitch/180;
    end

    if nargin < 8
        av_area = [10 10];
    end

    is2d = size(X,2)>1;
    if is2d
        X = X';
        Y = Y';
        Z = Z';
        nshapes = size(X,2);
    else
        nind = [0;find(isnan(X))];
        nshapes = length(nind)-1;
    end

    im = zeros(orig_im_size);

    [az1,el] = cart2sph(X-x,Y-y,Z-z);
    az = 1+mod((orig_im_size(2)-1)*(th-az1+pi)/(2*pi),orig_im_size(2)-1);
    el = orig_im_size(1)*(0.5+(pitch-el)/el_max);
    if is2d
        cazmid_all = 1+mod((orig_im_size(2)-1)* ...
            (th-atan2((Y+circshift(Y,-1))/2-y,(X+circshift(X,-1))/2-x)+pi)/(2*pi), ...
            orig_im_size(2)-1);
    else
        az = az(:)';
        el = el(:)';
    end

    for i = 1:nshapes
    %     fprintf('i: %d\n',i);
        if is2d
            cind = sub2ind(size(X),1:size(X,1),i*ones(1,size(X,1)));
            cazmid = cazmid_all(cind);
        else
            cind = nind(i)+1:nind(i+1)-1;
            cX = X(cind)-x;
            cY = Y(cind)-y;
            cazmid = mod(orig_im_size(2)* ...
                (th-atan2((cY+circshift(cY,-1))/2,(cX+circshift(cX,-1))/2)+pi)/(2*pi), ...
                orig_im_size(2));
        end

        caz = [az(cind),az(cind(1))];
        cel = [el(cind),el(cind(1))];

        rcel = round(el(cind));
        if all(rcel < 1) || all(rcel > orig_im_size(1))
            continue;
        end

        fullaz = [];
        fullel = [];
        for j = 1:length(caz)-1
            if caz(j)==caz(j+1)
                cazsign = 0;
            elseif (caz(j) < cazmid(j) && cazmid(j) < caz(j+1)) || ...
                   (cazmid(j) < caz(j+1) && caz(j+1) < caz(j)) || ...
                   (caz(j+1) < caz(j) && caz(j) < cazmid(j))
                cazsign = 1;
            else
                cazsign = -1;
            end

            if sign(caz(j)-caz(j+1))==cazsign
                caznxt = caz(j+1)+cazsign*(orig_im_size(2)-1);
            else
                caznxt = caz(j+1);
            end

            azdiff = caznxt-caz(j);
            eldiff = cel(j+1)-cel(j);
            if abs(azdiff) > abs(eldiff)
                if azdiff==0
                    fullaz2 = [caz(j),caz(j)];
                else
                    fullaz2 = fillvals(caz(j),caznxt,cazsign);
                end
                fullaz = [fullaz,fullaz2];

                fullel2 = linspace(cel(j),cel(j+1),length(fullaz2));
                fullel = [fullel,fullel2];

    %             fprintf('1\naz: %d\nel: %d\n\n',length(fullaz),length(fullel));
            else
                if eldiff==0
                    fullel2 = [cel(j),cel(j)];
                else
                    fullel2 = fillvals(cel(j),cel(j+1),sign(eldiff));
                end
                fullel = [fullel,fullel2];

                fullaz2 = linspace(caz(j),caznxt,length(fullel2));
                fullaz = [fullaz,fullaz2];

    %             fprintf('2\naz: %d\nel: %d\n\n',length(fullaz),length(fullel));
            end

    %         cim(sub2ind(size(cim),fullel2,1+mod(round(fullaz2),size(im,2)))) = true;
    %         figure(1);clf
    %         imshow(cim)
    %         disp(j)
    %         pause(1);
        end

    %     if any(5-fullel==2)
    %         keyboard
    %     end

    %     disp([range(fullel) range(fullaz)])
%         if all(fullel==fullel(1)) || all(fullaz==fullaz(1)) % zero width
%             continue;
%         end
        del = diff(fullel);
        daz = diff(fullaz);
        sel = del~=0 | daz~=0;
%         del = del(sel);
%         daz = daz(sel);

        fullel = 1+fullel(sel);
        fullaz = 1+mod(fullaz(sel),orig_im_size(2));
        [elflo,elfhi] = fllohi(fullel);
        [azflo,azfhi] = fllohi(fullaz);
        uarr = unique([elflo,elfhi; azflo, azfhi]','rows');
        lows = min(uarr);
%         uarr = uarr-lows+1;
        
        csz = [range(uarr(:,1)),range(uarr(:,2))];
        if any(csz==0)
            continue;
        end
        cim = zeros(csz);
        for j = 1:size(uarr,1)
            cuaz = fullaz(azflo==uarr(j,2) | azfhi==uarr(j,2));
            cuel = fullel(elflo==uarr(j,1) | elfhi==uarr(j,1));
            cim(uarr(j,1)-lows(1)+1,uarr(j,2)-lows(2)+1) = range(cuel)*range(cuaz);
%             keyboard
        end
        
%         loel = min(fullel);
%         hiel = max(fullel);
%         loaz = min(fullaz);
%         hiaz = max(fullaz);
% 
%     %     disp(max(orig_im_size(1)-fullel+1))
%         csz = [hiel-loel,hiaz-loaz+1];
%         if any(csz==0)
%             continue;
%         end
%         cim = false(csz); %+mod(-csz,1/fillthresh));
%         celind = fullel-loel;
%         cazind = fullaz-loaz+1;
%         sel = celind~=0;
%         cim(sub2ind(size(cim),celind(sel),cazind(sel))) = true;

    %     figure(1);clf
    %     image(cim)
    %     colormap gray
    %     pause(1);

        eli = lows(1)-1+(0:size(cim,1)-1);
%         eli(eli>size(im,1)) = [];
        azi = lows(2)+(0:size(cim,2)-1);
        azi(azi==0) = size(im,2);
        im(eli,azi) = min(1,im(eli,azi)+cim);

    %     figure(10);clf
    %     imshow(im)
    %     drawnow
    %     pause(1);
    end
    
    im = 1-im;

    horizel = round(size(im,1)*max(0,el_max/2 - pitch)/el_max);
    if horizel>0 || ~isempty(av_area)
        im = im2double(im);

        horizrows = (size(im,1)-horizel):size(im,1);
        im(horizrows,:) = max(0,im(horizrows,:)-0.5);

        if ~isempty(av_area)
            im=blockproc(im,av_area,@(x)mean2(x.data));
        end
    end
end

function vals=fillvals(from,to,sgn)
    vals = [];
    switch sgn
        case 1
            vals = ceil(from):floor(to);
        case -1
            vals = floor(from):-1:ceil(to);
    end
    if isempty(vals)
        vals = [from,to];
    else
        if vals(1)~=from
            vals = [from,vals];
        end
        if vals(end)~=to
            vals = [vals,to];
        end
    end
end

function [uarr,arr,ind]=flarr(vals)
    arr = floor(vals);
    rind = arr==vals;
    arr = [arr,arr(rind)-1];
    ind = [1:length(vals), find(rind)];
    uarr = unique(arr);
end

function [flo,fhi]=fllohi(vals)
    [flo,fhi] = deal(floor(vals));
    ind = flo==vals;
    flo(ind) = flo(ind)-1;
end