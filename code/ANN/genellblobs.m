% new on 20150505:
% COMBINED VERSIONS OF ELAZ AND ORSIEL NETS + NEW "MULTIAZ" NET (ORSIEL
% OVER AZES)

function genellblobs(whnet,apollo_starti)
% constants
debug = false;
netsuffixes = { 'elaz', 'orsiel+az', 'orsiel', 'elaz+or' };
netnvars = [ 2, 3, 3, 2 ];
if debug
    disp('DEBUG MODE, NOT SAVING')
end

% tweakables
pgeb.approxnvar = 1e4;
pgeb.nwave = 2;
pgeb.maxfreq = 30;
pgeb.maxamp = 1;
pgeb.scalemax = 1;
pgeb.a2b = 0.5;
pgeb.zrrng = [0.05 0.125];
pgeb.im_size = [120 270];
pgeb.blob_im_size = [120 360];
% pgeb.orsiel_az = -120:30:0; % for orsiel+az
pgeb.orsiel_az = linspace(-135,135,22);
pgeb.elaz_or = linspace(0,90,5);
pgeb.useblackblobs = true;
pgeb.usenothreshkerns = true;
pgeb.rim_size = [2 7; 4 7; 3 14];

rim_px_n = prod(pgeb.rim_size,2);
rim_whpx = cumsum([0; rim_px_n]);
xoff = (pgeb.blob_im_size(2)-pgeb.im_size(2))/2;

% apollo consts
if nargin < 2
    apollo_starti = 1;
end
pgeb.apollo_starti = apollo_starti;

pgeb.apollo_job = str2double(getenv('JOB_ID'));
if isnan(pgeb.apollo_job)
    pgeb.apollo_job = 0;
end

% load kernels; resize
if pgeb.usenothreshkerns
    load('vf_kernels_nothresh.mat','vf_avkernels*');
    kerns = [vf_avkernels_r2,vf_avkernels_r4];
    rkerns = resizekernel_nothresh(kerns,pgeb.im_size);
else
    load('vf_kernels.mat','vf_avkernels*');
    kerns = [vf_avkernels_r2,vf_avkernels_r4];
    rkerns = resizekernel(kerns,pgeb.im_size,0.25);
end

if ~debug
    dname = fullfile(mfiledir,'../../data/ANNs');
    if ~exist(dname,'dir')
        mkdir(dname)
    end
end
for i = 1:length(whnet)
    fprintf('----------------------------\n%d/%d: %s\n',i,length(whnet),netsuffixes{whnet(i)})
    
    pgeb.ndataeach = round(pgeb.approxnvar.^(1./netnvars(whnet(i))));
    ndata = pgeb.ndataeach.^netnvars(whnet(i));
    fprintf('ndata: %d (%d^%d)\n',ndata,pgeb.ndataeach,netnvars(whnet(i)));
    
    %% parameters
    el = linspace(-pgeb.im_size(1)/2,pgeb.im_size(1)/2,pgeb.ndataeach);
    if whnet(i)==1 || whnet(i)==4 % elaz
        areas = pi*pgeb.a2b*30.^2;
        az = linspace(-pgeb.im_size(2)/2,pgeb.im_size(2)/2,pgeb.ndataeach);
        if whnet(i)==4
            orients = pgeb.elaz_or;
        else
            orients = 0;
        end
    else
        areas = linspace(0,pi*pgeb.a2b*60.^2,pgeb.ndataeach+1);
        areas = areas(2:end);
        orients = linspace(0,90,pgeb.ndataeach);
        
        if whnet(i)==2 % orsiel+az
            az = pgeb.orsiel_az;
        else % orsiel
            lk = kerns(cell2mat({kerns.isleft}));
            cent = mean( cell2mat({lk.cent}') ./ ...
            cell2mat(cellfun(@(x)[size(x,2),size(x,1)],{lk.k},'UniformOutput',false)') );
            az = round(pgeb.im_size(2)*(cent(1)-0.5));
        end
    end
    
    [areas,orients,el,az] = ndgrid(areas,orients,el,az);
    areas = areas(:); orients = orients(:); el = el(:); az = az(:);
    ndata = length(areas);
    
    majoraxis = sqrt(areas./(pgeb.a2b*pi));
    minoraxis = majoraxis*pgeb.a2b;
    
    [truearea,aerr,scale] = deal(NaN(ndata,1));
    [amp,freq,phi] = deal(NaN(ndata,pgeb.nwave));
    x_im = NaN([ndata,rim_whpx(end)]);
    x_kern = NaN([ndata,length(kerns)]);
    if ~debug
        startprogbar(50,ndata);
    end
    c = 1;
    while c <= ndata
        while true
            scale(c) = pgeb.scalemax*rand;
            camp = rand(1,pgeb.nwave);
            amp(c,:) = range(pgeb.zrrng)*camp/sum(camp);
            freq(c,:) = randi(pgeb.maxfreq,1,pgeb.nwave);
            phi(c,:) = rand(1,pgeb.nwave)*2*pi;
            
            blob = bwtrim(~ellblob(scale(c),amp(c,:),freq(c,:),phi(c,:),majoraxis(c),minoraxis(c),pi*orients(c)/180,pgeb.blob_im_size));
            
            rp = regionprops(blob,'Orientation','Area','MajorAxisLength','MinorAxisLength');
            if length(rp)==1 % then one blob has been found (no split/too small blobs)
                aerr(c) = 1-areas(c)./(pi*rp.MajorAxisLength*rp.MinorAxisLength);
                if debug
                    fprintf('a1: %f; b1: %f; rat: %f\na2: %f; b2: %f; rat: %f\ndiff: %f\n\n', ...
                            majoraxis(i),minoraxis(i),minoraxis(i)/majoraxis(i), ...
                            rp.MajorAxisLength,rp.MinorAxisLength,rp.MinorAxisLength/rp.MajorAxisLength,aerr(i));
                end
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
        im = false(pgeb.blob_im_size);
        try
            ind = sub2ind(pgeb.blob_im_size,min(pgeb.blob_im_size(1),round(pgeb.blob_im_size(1)/2+ys)),round(pgeb.blob_im_size(2)/2+xs));
        catch
            %             disp('sub2ind')
            continue
        end
        im(ind) = true;
        im = imrotate(im,orients(c)-rp.Orientation,'crop');
        im = shiftmat(im,0,round(el(c)));
        im = circshift(im,[0,round(az(c))]);
        im = im(:,xoff+(1:pgeb.im_size(2)));
        
        truearea(c) = sum(im(:));
        
        if debug
            figure(1);clf
            imshow(im)
            keyboard
        end
        
        %         if debug
        % %             figure(1);clf
        %
        %             [xs,ys] = meshgrid(linspace(-im_size(2)/2,im_size(2)/2,im_size(2)),linspace(-im_size(1)/2,im_size(1)/2,im_size(1)));
        %             sn = sind(orients(i));
        %             cs = cosd(orients(i));
        %             cy = ys-el(i);
        %             ell = hypot((xs*cs-cy*sn)/majoraxis(i),(xs*sn+cy*cs)/minoraxis(i)) <= 1;
        %
        %             alsubplot(1,1)
        %             imshow(ell);
        %
        %             alsubplot(3,1)
        %             imshow(im);
        %
        %             keyboard
        %         end
        if pgeb.useblackblobs
            im = ~im;
        end
        
        x_kern(c,:) = getacts(im,rkerns);
        
        for k = 1:length(rim_px_n)
            rim = imresize(im,pgeb.rim_size(k,:),'bilinear');
            x_im(c,rim_whpx(k)+1:rim_whpx(k+1)) = rim(:);
        end
        
        if ~debug && progbar
            return;
        end
        
        c = c+1;
    end

    if ~debug
        fname = sprintf('%s/ellblob_%s_job%07d_ind%05d.mat',dname,netsuffixes{whnet(i)},pgeb.apollo_job,pgeb.apollo_starti);
        savemeta(fname,'x_im','x_kern','el','az','orients','areas','pgeb','majoraxis','minoraxis','aerr');
    end
end
