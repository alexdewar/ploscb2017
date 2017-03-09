function nsd_prbump(dosave)
    if ~nargin
        dosave = false;
    end
    
    dname = fullfile(mfiledir,'../../data/antoinestim/touse');
    fnames = { '09_4_00_162_1_-107_triangles_com.png','09_2_00_162_1_+067_triangles_hollow_com.png' };
    im_size = [120 360];
    fov = 270;
    offrad = 10;
    
    knum = [1,2,4];
    
    newsz = [im_size(1), im_size(2)*fov/360];

    load('vf_kernels','vf_avkernels_r2');
    kerns = vf_avkernels_r2;
    
    ks = cell(1,length(knum));
    
    ks{1} = resizekernel(kerns,newsz,.25);
    nkerns = size(ks{1},3)*knum;
    
    kcents = cell2mat({kerns.cent}');
    ksz = cell2mat(cellfun(@size,{kerns.k},'UniformOutput',false)');
    rcents = 1+bsxfun(@times,(kcents-1)./(ksz(:,[2 1])-1),newsz([2 1]));

    canplace = true(newsz);
    [yi,xi] = ndgrid(1:newsz(1),1:newsz(2));
    for j = 1:size(rcents,1)
        exrad(rcents(j,:));
    end
    
    function exrad(c)
        canplace = canplace & hypot(yi-c(2),xi-c(1)) > offrad;
    end

    canplaceorig = canplace;

%     figure(1);clf
%     imshow(canplace)
%     keyboard
    for i = 2:length(knum)
        canplace = canplaceorig;
        szfac = 1/sqrt(knum(i));
        
        ks{i} = zeros(size(ks{1},1),size(ks{1},2),nkerns(i));
        for j = 1:nkerns
            ck = kerns(j).k;
            ccent = kerns(j).cent;

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

            addexkern(sck,rcents(j,:),scent,j,i);

            for k = 1:knum(i)-1
                placei = find(canplace);
                [ky,kx] = ind2sub(newsz,placei(randi(length(placei))));
                addexkern(sck,[kx ky],scent,k*nkerns(1)+j,i);
                exrad([kx ky]);
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
    
    function addexkern(kern,newcent,oldcent,ind,whkern)
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
        
        ks{whkern}(:,:,ind) = ckern;
    end

    figure(1);clf
    alsubplot(1,length(fnames),1,1);
    for i = 1:length(fnames)
        cim = im2double(rgb2gray(imread(fullfile(dname,fnames{i}))));
        
%         alsubplot(1,i)
%         imshow(cim);
        
        rim = circshift(cim,[0 round(size(cim,2)*90/360)]);
%         alsubplot(2,i)
%         imshow(rim);
        
        xdiff = (im_size(2)-newsz(2))/2;
        
        tcim = cim(:,xdiff+1:end-xdiff);
        trim = rim(:,xdiff+1:end-xdiff);
        diffs = NaN(1,length(ks));
        for j = 1:length(ks)
            diffs(j) = getRMSdiff(getacts(tcim,ks{j}),getacts(trim,ks{j}));
        end

        alsubplot(1,i)
        bar(diffs,'w')
        set(gca,'XTickLabel',nkerns)
        xlabel('number of RFs')
        ylabel('RMS difference')
        ylim([0 .25])
        
        andy_setbox
    end
    
    if dosave
        savefig('nsd_prbump.pdf',[15, 5])
    end
end
