function genellblobs(starti,nvar,nwave,maxfreq,maxamp,job)
    if ~nargin
        debug = true;
        
        starti = 1;
        nvar = 1000;
        nwave = 2;
        maxfreq = 30;
        maxamp = 1;
        job = 1;
    else
        debug = false;
    end

    a2b = 0.5;

    scale = rand(nvar,1);
    amp = rand(nvar,nwave);
    phi = 2*pi*rand(nvar,nwave);
    freq = randi(maxfreq,nvar,nwave);
    % refim = ~ellblob(0,2,0,majoraxis,minoraxis,im_size);
    % refrp = regionprops(refim,'MajorAxisLength','MinorAxisLength');
    % refmja = refrp.MajorAxisLength;
    % refmna = refrp.MinorAxisLength;

    % [mja,mna] = deal(NaN(nvar,1));


% % % % % %     function genwaveblobs2(ind)
    ind = starti+(0:nvar-1);

    load('ml_inputs_azimuths','a1','rots1','chts1');
%     ind = wh*ngen+(1:ngen);

    area = a1(ind);
    orient = rots1(ind);
    el = chts1(ind);
    
    majoraxis = sqrt(area./(a2b*pi));
    minoraxis = majoraxis*a2b;
        
%     maxfreq = 7;
    zrrng = [0.05 0.5]; %[0.05 0.125];

    im_size = [120 270];
    av_area = 30;
%     nwave = 5;

    load('vf_tkernels.mat','tkerns_r2','tkerns_r4');
    tkerns = [tkerns_r2,tkerns_r4];

%     zerorad = NaN(nvar,1);
    % orient = 180*rand(ngen,1);
    [truearea,aerr] = deal(NaN(nvar,1));
    [amp,freq,phi] = deal(NaN(nvar,nwave));
    x_im = NaN([nvar,prod(im_size./av_area)]);
    x_kern = NaN([nvar,length(tkerns)]);
    if ~debug
        startprogbar(1,nvar);
    end
    for i = 1:nvar
        while true
%             zerorad(i) = rand*range(zrrng)+zrrng(1);
            camp = rand(1,nwave);
            amp(i,:) = range(zrrng)*camp/sum(camp);
            freq(i,:) = randi(maxfreq,1,nwave);
            phi(i,:) = rand(1,nwave)*2*pi;

            % amp(i,:)
            blob = bwtrim(ellblob(amp(i,:),1,freq(i,:),phi(i,:),majoraxis(i),minoraxis(i),pi*orient(i)/180,im_size));

            try
                rp = regionprops(blob,'Orientation','Area','MajorAxisLength','MinorAxisLength');
                aerr(i) = 1-area(i)./(pi*rp.MajorAxisLength*rp.MinorAxisLength);
            catch
                rp.Area
                continue;
            end
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
%             disp('sub2ind')
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
        im = imrotate(im,orient(i)-rp.Orientation,'crop');
        im = shiftmat(im,0,round(el(i)));

        truearea(i) = sum(im(:));

        if debug
%             figure(1);clf
            
            [xs,ys] = meshgrid(linspace(-im_size(2)/2,im_size(2)/2,im_size(2)),linspace(-im_size(1)/2,im_size(1)/2,im_size(1)));
            sn = sind(orient(i));
            cs = cosd(orient(i));
            cy = ys-el(i);
            ell = hypot((xs*cs-cy*sn)/majoraxis(i),(xs*sn+cy*cs)/minoraxis(i)) <= 1;
            
            alsubplot(1,1)
            imshow(ell);
            
            alsubplot(3,1)
            imshow(im);
            
            pause(0.5);
%             keyboard
        end

        [act,tkerns] = getneuronactivations(~im,tkerns);
        x_kern(i,:) = act';

        rim = blockproc(im,av_area*[1 1],@(x)mean2(x.data));

        x_im(i,:) = rim(:);

        if ~debug && progbar
            return;
        end
    end
% end
%     im = ~ellblob(scale(i),amp(i,:),freq(i,:),phi(i,:),majoraxis,minoraxis,thoff(i),im_size);
% 
if ~debug
    fname = sprintf('%s/ellblob%05d_%05d.mat',mfiledir,job,starti);
    fprintf('Saving to %s...\n',fname);
    save(fname,'x_im','x_kern','majoraxis','minoraxis','im_size','scale','amp','phi','freq','maxfreq','maxamp','nvar','nwave','job','starti');
end