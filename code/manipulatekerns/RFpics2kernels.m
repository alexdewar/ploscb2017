 function RFpics2kernels(dosave)
 if ~nargin
     dosave = false;
 end
%     close all
 
    crng = [1 -1];
    thresh = 0.25;
    avthresh = 0.25;
    graybarthresh = 0.03; % std <=
    r2dir = [mfiledir '/../../data/receptive_fields_pics'];
    r4dir = [r2dir '/r4d'];
    toosmall = 20;
    fn = [mfiledir '/../../data/all/vf_kernels_new.mat'];

%     areas = [];
    cbim = mean(im2double(imread('colorbar.ppm')),2);
%     cbim = cbim(end:-1:1,:,:);
    cvals = linspace(crng(1),crng(2),size(cbim,1));
%     cbim(end+1,1,:) = [189.6717171717172,189.3888888888889,194.5101010101010];
%     cvals(end+1) = 0;
    [ikr4,imsr4,fn4,gn4] = getikernels(r4dir);
    [vf_kernels_r4,vf_avkernels_r4] = centerkernel(ikr4,imsr4,fn4,gn4);
    neuroncolormap = flipud(cbim(:,:));
%     pause
%     for knum=1:7
%         figure(1);clf
%         imagesc(avkernels_r4{knum})
%         colormap(neuroncolormap);
%         pause(1);
%     end
%     return
    [ikr2,imsr2,fn2,gn2] = getikernels(r2dir);
    [vf_kernels_r2,vf_avkernels_r2] = centerkernel(ikr2,imsr2,fn2,gn2);
%     pause
%     for knum=1:length(avkernels_r6)
%         figure(1);clf
%         imagesc(sign(avkernels_r6{knum}))
%         colormap(neuroncolormap);
%         pause(1);
%     end
%     keyboard
%     imsize_r4 = [size(imsr4,1) size(imsr4,2)];
%     imsize_r2 = [size(imsr2,1) size(imsr2,2)];
    if dosave
        disp('Saving...')
        save(fn,'vf_kernels_r4','vf_avkernels_r4','vf_kernels_r2', ...
            'vf_avkernels_r2','neuroncolormap', ...
            'toosmall','thresh','avthresh','graybarthresh');
    end

    function [kstr,avkstr]=centerkernel(ks,ims,flynum,glomnum)
        kstr(length(flynum)*2) = struct('k',[],'glomnum',[],'flynum',[], ...
             'cent',[],'isleft',[]);
%         avcents = [];
        for j = 1:length(flynum)
            kstr(j).glomnum = glomnum(j);
            kstr(j).flynum = flynum(j);
            kstr(j).isleft = true;
            j2 = j+length(flynum);
            kstr(j2).glomnum = glomnum(j);
            kstr(j2).flynum = flynum(j);
            kstr(j2).isleft = false;
            [kstr(j).k,kstr(j).cent,kstr(j2).k,kstr(j2).cent] ...
                = finishkernel(ks(:,:,j));
            
%             if glomnum(j)==5
% %                 figure(j);clf
% %                 showkernel(kstr(j));
%                 avcents(end+1,:) = kstr(j).cent;
%             end
        end

        gns = unique(glomnum)';
        avkstr(length(gns)*2) = struct('k',[],'glomnum',[],'flynum',[], ...
                        'cent',[],'isleft',[]);
        klefts = kstr(cell2mat({kstr.isleft}));
        for j = 1:length(gns)
            cglom = find(cell2mat({klefts.glomnum})==gns(j));
            avcent = mean(cell2mat({klefts(cglom).cent}'));
            
            fprintf('g%02d: %f,%f\n',gns(j),avcent);
            
%             if gns(j)==5
%                 keyboard
%             end
            asz = [size(ims,1) size(ims,2)];
            avk = NaN([asz,length(cglom)]);

%             figure(1);clf;hold on
            for k = 1:length(cglom)
                tk = ims(:,:,cglom(k));

                ccent = klefts(cglom(k)).cent;
                tk = shiftmat(tk,avcent(1)-ccent(1),avcent(2)-ccent(2));
                
                if gns(j)==11
                    figure(1)
                    subplot(2,3,k)
                    showkernel(sign(tk).*(abs(tk)>=thresh));
                    hold on
                    plot(avcent(1),avcent(2),'k+')
                end
%                 keyboard
                
                avk(:,:,k) = changematsize(tk,asz);
                
                
%                 figure(1)
%                 currentcentre = biggestarea(avk(:,:,k)>.25);
%                 plot(currentcentre(1),currentcentre(2),'g+');
%                 [curx,cury] = bw2polygon(avk(:,:,k));
%                 alfill(curx,cury,'r','facealpha',1/length(cglom));
%                 keyboard
            end
            
            cavk = nanmean(avk,3);
            cavk = sign(cavk).*(normalizerange(abs(cavk))>=avthresh);
            
            fns = cell2mat({klefts(cglom).flynum});
            
            avkstr(j).glomnum = gns(j);
            avkstr(j).flynum = fns;
            avkstr(j).isleft = true;
%             avkstr(j).offset = [NaN NaN];
            j2 = j+length(gns);
            avkstr(j2).glomnum = gns(j);
            avkstr(j2).flynum = fns;
            avkstr(j2).isleft = false;
%             avkstr(j2).offset = [NaN NaN];
            [avkstr(j).k,avkstr(j).cent,avkstr(j2).k,avkstr(j2).cent] ...
                = finishkernel(cavk,[size(ims,1) size(ims,2)]);
        end
    end

    function [k1,cent1,k2,cent2]=finishkernel(x,sz)        
        k1 = imfill(x==1,'holes')-imfill(x==-1,'holes');
        k1 = removetoosmall(k1);
        
        if nargin==2
            k1 = changematsize(bwtrim(k1),sz);
            
%             ksz = size(k1);
%             k1 = padarray(k1,floor((sz-ksz)/2));
%             k1 = padarray(k1,sz-size(k1),'post');
        else
            sz = [size(k1,1) size(k1,2)];
        end
%         off1 = [xtrim(1) ytrim(1)];
%         off2 = [xtrim(2) ytrim(2)];
        
        cent1 = biggestarea(k1==1);
        cent2 = [sz(2)-cent1(1) cent1(2)];
        
        inh = k1==-1;
        k1(inh) = -1/sum(inh(:));
        exc = k1==1;
        k1(exc) = 1/sum(exc(:));
        k2 = fliplr(k1);
    end

    function x=removetoosmall(x)
        cbw = bwlabeln(x==-1);
        cexc = bwlabeln(x==1);
        cbw(cexc>0) = max(cbw(:))+cexc(cexc>0);
        rp = regionprops(cbw,'Area');
        for k = 1:length(rp)
%             areas(end+1) = rp(k).Area;
            if rp(k).Area <= toosmall
                x(cbw==k) = 0;
            end
        end
    end

    function cent=biggestarea(x)
        stats = regionprops(bwlabeln(x),'Area','Centroid');
        maxarea = 0;
        whmax = NaN;
        for j = 1:length(stats)
            if stats(j).Area > maxarea
                maxarea = stats(j).Area;
%                 areas(end+1) = stats(j).Area;
                whmax = j;
            end
        end
        cent = stats(whmax).Centroid;
    end

    function [ks,ims,flynum,glomnum] = getikernels(dname)
        owd = pwd;
        cd(dname);
        ikfn = 'ikernels.mat';
        if exist(ikfn,'file')
            load(ikfn);
        else
            rd = dir('*.jpg');
            flynum = NaN(length(rd),1);
            glomnum = NaN(length(rd),1);
            for i = 1:length(rd)
                name = rd(i).name;
                disp(name)
                glomnum(i) = str2double(name(2:3));
                flynum(i) = str2double(name(5));
                im = imread(name);
                if ~exist('ims','var')
                    sz = size(im);
%                     radperpix = fov./sz(1:2);
                    [ims,ks] = deal(NaN([sz(1:2) length(rd)]));
                end
                im = im2double(im);
                gbt = std(mean(im,3),[],2)<=graybarthresh;
%                 if any(gbt)
%                     figure(1);clf
%                     tim = im;
%                     subplot(1,2,1);
%                     imagesc(tim);
%                     subplot(1,2,2);
%                     tim(gbt,:,:) = 0;
%                     imagesc(tim);
%                     keyboard
%                 end
                im(gbt,:,:) = 1;
                cval = im2cval(im);
                ims(:,:,i) = cval;
                ks(:,:,i) = cval2kern(cval);
            end
            save(ikfn,'ks','ims','flynum','glomnum','radperpix');
        end
    %     areas = [];
        cd(owd);
    end

    function im2=im2cval(im)
%         im = im2double(im);
        im2 = NaN(size(im,1),size(im,2));
        for j = 1:size(im,1)
            for k = 1:size(im,2)
                mdiff = mean(abs(bsxfun(@minus,im(j,k,:),cbim)),3);
                [mval,mind] = min(mdiff);
                im2(j,k) = cvals(mind);
            end
        end
    end

    function k=cval2kern(cval)
        k = sign(cval).*(abs(cval) >= thresh);
    end
end