function im=getviewfast(x,y,z,th,X,Y,Z,av_area,el_max,orig_im_size,pitch)
% function im=getviewfast(x,y,z,th,X,Y,Z,av_area,el_max,orig_im_size,pitch)
%     minshapesize = 1;
    
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
        pitch = pitch*(pi/180);
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

    im = false(orig_im_size);

    [az1,el] = cart2sph(X-x,Y-y,Z-z);
    az = mod((orig_im_size(2)-1)*(th-az1+pi)/(2*pi),orig_im_size(2)-1);
    el = (orig_im_size(1)-1)*(0.5+(pitch-el)/el_max);
%     if ~is2d
%         cazmid_all = 1+mod((orig_im_size(2)-1)* ...
%             (th-atan2((Y+circshift(Y,-1))/2-y,(X+circshift(X,-1))/2-x)+pi)/(2*pi), ...
%             orig_im_size(2)-1);
%     else
        az = az(:)';
        el = el(:)';
%     end

    for i = 1:nshapes
    %     fprintf('i: %d\n',i);
        if is2d
            cind = sub2ind(size(X),1:size(X,1),i*ones(1,size(X,1)));
%             cazmid = cazmid_all(cind);
        else
            cind = nind(i)+1:nind(i+1)-1;
%             cX = X(cind)-x;
%             cY = Y(cind)-y;
%             cazmid = mod(orig_im_size(2)* ...
%                 (th-atan2((cY+circshift(cY,-1))/2,(cX+circshift(cX,-1))/2)+pi)/(2*pi), ...
%                 orig_im_size(2));
        end
        
        cel = round(el(cind));
        his = cel > orig_im_size(1)-1;
        los = cel < 0;
%         if all(hiorlo)
%             continue;
%         end
        
        hi2hi = his & his([2:end,1]);
        lo2lo = los & los([2:end,1]);
        elouts = ~hi2hi & ~lo2lo;
        cel = cel(elouts);
        cind = cind(elouts);
        caz = round(az(cind));
        
        cel = max(0,min(cel,orig_im_size(1)-1));
        
%         caz1 = az1(cind);
%         caz2 = caz1([2:end,2]);
        cazsign = sign(caz([2:end,1])-caz); %sign(caz2-caz1);
        
%         cazsel = caz1~=0 & sign(caz1) == -sign(caz2);
%         if any(cazsel)
%             y1 = Y(cind);
%             y2 = [y1(2:end);y1(1)];
%             cazsign(cazsel) = sign(y2(cazsel)-y1(cazsel));
%         end

        celsign = sign(cel([2:end,1])-cel);
        
        cim = false(orig_im_size); %false(range(cel)+1,range(caz)+1);
        del = cel([2:end,1])-cel;
        derr = abs(del./(caz([2:end,1])-caz));
        caz = 1+[caz,caz(1)];
        cel = 1+[cel,cel(1)];
        for j = 1:length(caz)-1
            if del(j)==0 % horizontal line or single pixel
                if cazsign(j)==0
                    cim(cel(j),caz(j)) = true;
                else
                    cim(cel(j),fillvals(caz(j),caz(j+1),cazsign(j),orig_im_size(2))) = true;
                end
            elseif isinf(derr(j)) % vertical line
                cim(fillvals(cel(j),cel(j+1),celsign(j),orig_im_size(1)),caz(j)) = true;
            else
                err = 0;
                y = cel(j);
                for x = fillvals(caz(j),caz(j+1),cazsign(j),orig_im_size(2))
                    cim(y,x) = true;
                    err = err+derr(j);
                    while err >= 0.5
                        cim(y,x) = true;
                        y = min(orig_im_size(1),y+sign(del(j)));
                        err = err-1;
                    end
                end
            end
            
%             if ~all(size(cim)==orig_im_size)
%                 keyboard
%             end
        end
        
        cimfill = false(size(cim));
        
        dcim = diff([false(1,orig_im_size(2));cim])==1;
        only1edge = sum(dcim)==1;
        
        cimfill(:,only1edge) = cim(:,only1edge);
        cimfill(:,~only1edge) = cim(:,~only1edge) | mod(cumsum(dcim(:,~only1edge)),2)==1;
        
%         figure(1);clf
%         subplot(2,1,1)
%         imshow(cimfill)
%         subplot(2,1,2)
%         imagesc(diff([false(1,orig_im_size(2));cim]))
%         colorbar
%         keyboard
        
        im = im | cimfill;
        
    end
    
    im = ~im;

    horizel = round(size(im,1)*max(0,el_max/2 - pitch)/el_max);
    if horizel>0 || ~isempty(av_area)
        im = im2double(im);

%         horizrows = (size(im,1)-horizel):size(im,1);
%         im(horizrows,:) = max(0,im(horizrows,:)-0.5);

        if ~isempty(av_area)
            im=blockproc(im,av_area,@(x)mean2(x.data));
        end
    end
end

function vals=fillvals(from,to,sgn,lim)
    if sgn==0
        vals = from;
    else
        truesgn = sign(to-from);
        if truesgn==sgn
            vals = from:sgn:to;
        elseif sgn==-1
            vals = [1:from, to:lim];
        else
            vals = [1:to, from:lim];
        end
    end
end
