 function RFpics2kernels_nothresh(dosave)
    if ~nargin
        dosave = false;
    end
%     close all
 
    crng = [1 -1];
    stdthresh = 0.5;
    avthresh = 0.25;
    graybarthresh = 0.03; % std <=
    r2dir = [mfiledir '/../../data/receptive_fields_pics'];
    r4dir = [r2dir '/r4d'];
    toosmall = 20;
    fn = [mfiledir '/../../data/all/vf_kernels_nothresh.mat'];

%     areas = [];
    cbim = mean(im2double(imread('colorbar.ppm')),2);
    neuroncolormap = flipud(cbim(:,:));
%     cbim = cbim(end:-1:1,:,:);
    cvals = linspace(crng(1),crng(2),size(cbim,1));
%     cbim(end+1,1,:) = [189.6717171717172,189.3888888888889,194.5101010101010];
%     cvals(end+1) = 0;
    [imsr4,fn4,gn4] = getikernels(r4dir,'r4');
    [vf_kernels_r4,vf_avkernels_r4] = centerkernel(imsr4,fn4,gn4);

%     pause
%     for knum=1:7
%         figure(1);clf
%         imagesc(avkernels_r4{knum})
%         colormap(neuroncolormap);
%         pause(1);
%     end
%     return
    [imsr2,fn2,gn2] = getikernels(r2dir,'r2');
    [vf_kernels_r2,vf_avkernels_r2] = centerkernel(imsr2,fn2,gn2);
%     pause
%     for knum=1:length(avkernels_r6)
%         figure(1);clf
%         imagesc(sign(avkernels_r6{knum}))
%         colormap(neuroncolormap);
%         pause(1);
%     end
%     keyboard
%     imsize_r4 = [size(imsr4,1) size(imsr4,2)];
% %     imsize_r2 = [size(imsr2,1) size(imsr2,2)];
    if dosave
        save(fn,'vf_kernels_r4','vf_avkernels_r4','vf_kernels_r2', ...
            'vf_avkernels_r2','neuroncolormap', ...
            'toosmall','stdthresh','avthresh','graybarthresh');
    else
        dump2base(true);
    end

    function [kstr,avkstr]=centerkernel(ims,flynum,glomnum)
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
            [kstr(j).k,kstr(j).cent,kstr(j2).k,kstr(j2).cent] ... % switched from ks to ims
                = finishkernel(ims(:,:,j));
            
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
        
%         figure(1);clf
        
        for j = 1:length(gns)
            cglom = find(cell2mat({klefts.glomnum})==gns(j));
            avcent = mean(cell2mat({klefts(cglom).cent}'));
            
            fprintf('g%02d: %f,%f\n',gns(j),avcent);
            
%             if gns(j)==5
%                 keyboard
%             end
            asz = [size(ims,1) size(ims,2)];
            avk = NaN([asz,length(cglom)]);

%             a=ims(:,:,cglom);
%             figure(1);clf
            for k = 1:length(cglom)
                tk = ims(:,:,cglom(k));
                ccent = klefts(cglom(k)).cent;
                
%                 subplot(2,1,1)
%                 imagesc(tk)
%                 colormap(neuroncolormap)
%                 hold on
%                 plot(ccent(1),ccent(2),'g+',avcent(1),avcent(2),'b+')

                tk = changematsize(shiftmat(tk,avcent(1)-ccent(1),avcent(2)-ccent(2)),asz);

                if gns(j)==11
                    figure(101)
                    subplot(2,3,k)
                    showkernel_nothresh(tk);
                    hold on
                    plot(avcent(1),avcent(2),'k+')
                end
                
%                 subplot(2,1,2)
%                 imagesc(tk)
%                 colormap(neuroncolormap)
%                 hold on
%                 plot(ccent(1),ccent(2),'g+',avcent(1),avcent(2),'b+')
%                 
%                 keyboard

                avk(:,:,k) = tk;
            end
            
            cavk = nanmean(avk,3);
%             cavk = sign(cavk).*(normalizerange(abs(cavk))>=avthresh);
            
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
            
%             figure(2);clf
%             imagesc(avkstr(j).k)
%             colormap(neuroncolormap)
%             keyboard
        end
    end

    function [k1,cent1,k2,cent2]=finishkernel(x,sz)        
%         k1 = imfill(x==1,'holes')-imfill(x==-1,'holes');
        
        k1 = removetoosmall(x);
        
        inh = k1<0;
        k1(inh) = -k1(inh)./min(k1(:));
        exc = k1>0;
        k1(exc) = k1(exc)./max(k1(:));
        
        if nargin==2
            k1 = changematsize(bwtrim(k1),sz);
            
%             figure(1);clf
%             subplot(2,1,1)
%             imagesc(k1)
%             colormap(neuroncolormap)
%             subplot(2,1,2)
%             imagesc(cval2kern(k1))
%             colormap(neuroncolormap)
%             keyboard
        else
            sz = [size(k1,1) size(k1,2)];
        end
        
        k1 = k1.*(abs(k1) > stdthresh*std(k1(:)));
        
        cent1 = biggestarea(k1>0);
        cent2 = [sz(2)-cent1(1) cent1(2)];
        
%         inh = k1==-1;
%         k1(inh) = -1/sum(inh(:));
%         exc = k1==1;
%         k1(exc) = 1/sum(exc(:));
        k2 = fliplr(k1);
    end

    function x=removetoosmall(x)
        cbw = bwlabeln(x<0);
        cexc = bwlabeln(x>0);
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

    function [ims,flynum,glomnum] = getikernels(dname,suffix)
        owd = pwd;
        cd(dname);
        ikfn = sprintf('%s/../../data/all/ikernels_%s_nothresh.mat',mfiledir,suffix);
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
                    ims = NaN([sz(1:2) length(rd)]);
                end
                im = im2double(im);
                gbt = std(mean(im,3),[],2)<=graybarthresh;

                im(gbt,:,:) = 1;
                cval = im2cval(im);
                ims(:,:,i) = cval;
            end
            save(ikfn,'ims','flynum','glomnum');
        end
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