function im=getviewfast(x,y,z,th,X,Y,Z,av_area,el_max,orig_im_size,pitch)
% function im=getviewfast(x,y,z,th,X,Y,Z,av_area,el_max,orig_im_size,pitch)
fillthresh = 0.2;

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

bigsz = orig_im_size/fillthresh;
im = false(bigsz);

[az1,el] = cart2sph(X-x,Y-y,Z-z);
az = 1+mod(bigsz(2)*(th-az1+pi)/(2*pi),bigsz(2));
el = bigsz(1)*(0.5+(pitch-el)/el_max);
if is2d
    cazmid_all = mod(bigsz(2)* ...
        (th-atan2((Y+circshift(Y,-1))/2-y,(X+circshift(X,-1))/2-x)+pi)/(2*pi), ...
        bigsz(2));
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
        cazmid = mod(bigsz(2)* ...
            (th-atan2((cY+circshift(cY,-1))/2,(cX+circshift(cX,-1))/2)+pi)/(2*pi), ...
            bigsz(2));
    end

    caz = [az(cind),az(cind(1))];
    cel = [el(cind),el(cind(1))];

    rcel = round(el(cind));
    if all(rcel < 1) || all(rcel > bigsz(1))
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
            caznxt = caz(j+1)+cazsign*(bigsz(2)-1);
        else
            caznxt = caz(j+1);
        end
        
        azdiff = caznxt-caz(j);
        eldiff = cel(j+1)-cel(j);
        if abs(azdiff) > abs(eldiff)
            if azdiff==0
                fullaz2 = [caz(j),caz(j)];
            else
                fullaz2 = caz(j):cazsign:caznxt;
                if fullaz2(end)~=caznxt
                    fullaz2(end+1) = caznxt;
                end
            end
            fullaz = [fullaz,fullaz2];

            fullel2 = linspace(cel(j),cel(j+1),length(fullaz2));
            fullel = [fullel,fullel2];
                
%             fprintf('1\naz: %d\nel: %d\n\n',length(fullaz),length(fullel));
        else
            if eldiff==0
                fullel2 = [cel(j),cel(j)];
            else
                fullel2 = cel(j):sign(eldiff):cel(j+1);
                if fullel2(end)~=cel(j+1)
                    fullel2(end+1) = cel(j+1);
                end
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

    fullel = round(fullel);
    fullaz = 1+mod(round(fullaz),bigsz(2));
    loel = min(fullel);
    hiel = max(fullel);
    loaz = min(fullaz);
    hiaz = max(fullaz);
    
%     disp(max(orig_im_size(1)-fullel+1))
    csz = [hiel-loel,hiaz-loaz+1];
    if any(csz==0)
        continue;
    end
    cim = false(csz+mod(-csz,1/fillthresh));
    celind = fullel-loel;
    cazind = fullaz-loaz+1;
    sel = celind~=0;
    cim(sub2ind(size(cim),celind(sel),cazind(sel))) = true;

%     figure(1);clf
%     image(cim)
%     colormap gray
%     pause(1);
    
    if size(cim,1)>1
        ecim = cim;
        ecim(:,sum(ecim)<2) = false;
        df = [diff([false(1,size(ecim,2));ecim(1:end-1,:)])==1;ecim(end,:)];
        df(:,sum(df)<2) = false;
        edges = cumsum(df);
        cim(mod(edges,2)==1) = true;
    end
    
%     disp(2-loel);
    cst = max(1,2-loel);
    cend = min(hiel,bigsz(1))-loel;
    cim = cim(cst:cend,:);
    
%     figure(1);clf
%     image(cim*255)
%     colormap gray
%     drawnow
%     pause(1);

%     keyboard

    eli = max(loel,1)+(0:size(cim,1)-1);
    azi = loaz+(0:size(cim,2)-1);
    im(eli,azi) = im(eli,azi) | cim; %blockproc(cim,[1 1]/fillthresh,@(x)any(x.data(:)));
    
%     figure(10);clf
%     imshow(im)
%     drawnow
%     pause(1);
end

im = ~im;

horizel = round(size(im,1)*max(0,el_max/2 - pitch)/el_max);
if horizel>0 || ~isempty(av_area)
    im = im2double(im);

    horizrows = (size(im,1)-horizel):size(im,1);
    im(horizrows,:) = max(0,im(horizrows,:)-0.5);

    if ~isempty(av_area)
        im=blockproc(im,av_area,@(x)mean2(x.data));
    end
end