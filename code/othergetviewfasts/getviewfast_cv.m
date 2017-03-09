function im=getviewfast_cv(x,y,z,th,X,Y,Z,av_area,el_max,orig_im_size,pitch)
% function im=getviewfast(x,y,z,th,X,Y,Z,av_area,el_max,orig_im_size,pitch)
%     fillthresh = 0.01;

    if nargin < 10 || isempty(orig_im_size)
        orig_im_size = [170 900];
    end
    if nargin < 9 || isempty(el_max)
        el_max = 5*pi/12;
    else
        el_max = pi*el_max/180;
    end

    if nargin < 11 || isempty(pitch)
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

    im = zeros(orig_im_size,'uint8');

    [az1,el] = cart2sph(X-x,Y-y,Z-z);
    az = 1+mod((orig_im_size(2)-1)*(th-az1+pi)/(2*pi),orig_im_size(2));
    el = 1+(orig_im_size(1)-1)*(0.5+(pitch-el)/el_max);
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
%             cazmid = cazmid_all(cind);
        else
            cind = nind(i)+1:nind(i+1)-1;
%             cX = X(cind)-x;
%             cY = Y(cind)-y;
%             cazmid = mod(orig_im_size(2)* ...
%                 (th-atan2((cY+circshift(cY,-1))/2,(cX+circshift(cX,-1))/2)+pi)/(2*pi), ...
%                 orig_im_size(2));
        end

        caz = az(cind);
        azlo = floor(min(caz));
        cel = el(cind);
        ello = floor(min(cel));
        ccoord = [1+caz-azlo;1+cel-ello];

%         rcel = round(cel);
%         if all(rcel < 1) || all(rcel > orig_im_size(1))
%             continue;
%         end
        
        cim = insertShape(zeros(ceil(max(caz))-azlo,ceil(max(cel))-ello,'uint8'), ...
                          'FilledPolygon',ccoord(:)','Color','white','Opacity',1);
        eli = ello+(0:size(cim,1)-1);
        azi = azlo+(0:size(cim,2)-1);
        try
            im(eli,azi) = im(eli,azi) + cim(:,:,1);
        catch
            disp(max(eli))
        end
%         keyboard
%         figure(1);clf
%         imshow(im)
%         drawnow
%         pause(.3)
    end
    
    im = 255-im;
    
%     figure(1);clf
%     imshow(im)
%     return

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