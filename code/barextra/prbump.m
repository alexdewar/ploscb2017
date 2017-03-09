function prbump
    nsim = 100;
    
    dname = './dispatchrecap/antoinestim/new/';
    fnames = { 'triangles.png','triangles_com.png','triangles_hollow.png','triangles_hollow_com.png' };
    im_size = [120 360];
    fov = 270;
    offrad = 10;
    
    knum = [1,1,2,4];
    
    newsz = [im_size(1), im_size(2)*fov/360];

    load('vf_kernels','vf_avkernels_r2');
    kerns = vf_avkernels_r2;
    
    ks = cell(1,length(knum));
    
    ks{1} = resizekernel(kerns,newsz,.25);
    nkerns = size(ks{1},3)*knum;
    
    kcents = cell2mat({kerns.cent}');
    ksz = cell2mat(cellfun(@size,{kerns.k},'UniformOutput',false)');
    rcents = 1+bsxfun(@times,(kcents-1)./(ksz(:,[2 1])-1),newsz([2 1]));

%     canplace = true(newsz);
    [yi,xi] = ndgrid(1:newsz(1),1:newsz(2));
%     for k = 1:size(rcents,1)
%         exrad(rcents(k,:));
%     end
    
    function exrad(c)
        canplace = canplace & hypot(yi-c(2),xi-c(1)) > offrad;
    end

%     canplaceorig = canplace;

%     figure(1);clf
%     imshow(canplace)
%     keyboard

    startprogbar(50,nsim.*sum(nkerns(2:end)),'generating kernels - ')
    for i = 2:length(knum)
        szfac = 1/sqrt(knum(i));
        
        ks{i} = zeros(size(ks{1},1),size(ks{1},2),nkerns(i),nsim);
        for j = 1:nsim
            canplace = true(newsz);
            
            for k = 1:nkerns(1)
                ck = kerns(k).k;
                ccent = kerns(k).cent;

                anyy = any(ck,2);
                anyx = any(ck);
                ybeg = find(anyy,1);
                yend = find(anyy,1,'last');
                xbeg = find(anyx,1);
                xend = find(anyx,1,'last');
                sck = ck(ybeg:yend,xbeg:xend);

                sck = resizekernel(sck,szfac,0.25);
                scent = szfac*[ccent(1)-xbeg+1,ccent(2)-ybeg+1];

        %         figure(2);clf
        %         showkernel(sck,scent)
        %         keyboard

    %             addexkern(sck,rcents(k,:),scent,k,i);

                for l = 0:knum(i)-1
                    placei = find(canplace);
                    [ky,kx] = ind2sub(newsz,placei(randi(length(placei))));
                    addexkern(sck,[kx ky],scent,l*nkerns(1)+k,j,i);
                    exrad([kx ky]);
                    
                    if progbar
                        return;
                    end
                end

        %         figure(1);clf
        %         subplot(1,4,1)
        %         showkernel(exkerns(:,:,i),rcents(i,:));
        %         subplot(1,4,2)
        %         showkernel(ks{1}(:,:,i),rcents(i,:));
        %         title(sprintf('%.2f ',rcents(i,:)))
        %         subplot(1,4,3)
        %         showkernel(kerns(i))
        %         title(sprintf('%.2f ',kerns(i).cent));
        %         subplot(1,4,4)
        %         showkernel(exkerns(:,:,i+nkerns),[kx ky]);
        %         keyboard
            end
        end
    end
    
    function addexkern(kern,newcent,oldcent,whkern,whsim,whset)
        npre = floor(newcent-oldcent);
        
        ckern = zeros(newsz);
        ckern(1:size(kern,1),1:size(kern,2)) = kern;
        ckern = shiftmat(ckern,npre(1),npre(2));
        
%         yind = max(1,npre(2)+1):min(newsz(1),npre(2)+ssz(1));
%         xind = max(1,npre(1)+1):min(newsz(2),npre(1)+ssz(2));
%         ckern = kern(1:length(yind),1:length(xind));
        pos = ckern>0;
        neg = ckern<0;
        ckern(pos) = 1/sum(pos(:));
        ckern(neg) = -1/sum(neg(:));
        
        ks{whset}(:,:,whkern,whsim) = ckern;
    end

    figure(1);clf
    alsubplot(3,length(fnames),1,1);
    startprogbar(50,length(fnames).*sum(cellfun(@(x)size(x,4),ks)),'getting activations - ');
    for i = 1:length(fnames)
        cim = im2double(rgb2gray(imread([dname,fnames{i}])));
        
        alsubplot(1,i)
        imshow(cim);
        
        rim = circshift(cim,[0 round(size(cim,2)*90/360)]);
        alsubplot(2,i)
        imshow(rim);
        
        xdiff = (im_size(2)-newsz(2))/2;
        
        tcim = cim(:,xdiff+1:end-xdiff);
        trim = rim(:,xdiff+1:end-xdiff);
        [diffs,errs] = deal(zeros(1,length(ks)));
        for j = 1:length(ks)
            cdiffs = NaN(1,size(ks{j},4));
            for k = 1:length(cdiffs)
                cdiffs(k) = getRMSdiff(getacts(tcim,ks{j}),getacts(trim,ks{j}));
                
                if progbar
                    return;
                end
            end
            diffs(j) = mean(cdiffs);
            errs(j) = stderr(cdiffs);
        end

        alsubplot(3,i)
        barerr(diffs,errs,'w')
        set(gca,'XTickLabel',{[]; num2str(nkerns'); []})
        xlabel('number of RFs')
        ylabel('RMS difference')
        ylim([0 .25])
    end
    dump2base
end
