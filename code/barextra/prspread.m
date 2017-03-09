clear

nrfrow = 2;

nfigrow = 5;

dname = './dispatchrecap/antoinestim/new/';
d = dir([dname '*.png']);
fnames = {d.name};
% fnames = { 'triangles.png','triangles_com.png','triangles_hollow.png','triangles_hollow_com.png' };
im_size = [120 360];
fov = 270;
offrad = 10;

knum = [1,2,4];

newsz = [im_size(1), im_size(2)*fov/360];

load('vf_kernels','vf_avkernels_r2');
kerns = vf_avkernels_r2;

ks = resizekernel(kerns,newsz,.25);
ks(:,:,:,2) = NaN;

kcents = cell2mat({kerns.cent}');
ksz = cell2mat(cellfun(@size,{kerns.k},'UniformOutput',false)');
cents = NaN([size(kcents), 2]);
cents(:,:,1) = 1+bsxfun(@times,(kcents-1)./(ksz(:,[2 1])-1),newsz([2 1]));

kyd = newsz(1)/(1+nrfrow);
kxd = newsz(2)/(1+length(kerns)/nrfrow);
kys = kyd:kyd:newsz(1)-kyd;
kxs = kxd:kxd:newsz(2)-kxd;
[mx,my] = meshgrid(kxs,kys);
mx = mx(:)';
my = my(:)';

[cinds,rinds] = deal(1:length(mx));

xd = bsxfun(@minus,mx,cents(:,1,1));
yd = bsxfun(@minus,my,cents(:,2,1));
hs = hypot(yd,xd);
[rows,cols] = deal(true(1,length(mx)));
for i = 1:length(mx)
    [mins,inds] = min(hs(rows,cols));
    [minv,mini] = min(mins);
    
    cind = cinds(mini);
    cols(cind) = false;
    cinds(mini) = [];
    
    rind = rinds(inds(mini));
    rows(rind) = false;
    rinds(inds(mini)) = [];
    
    cx = mx(cind);
    cy = my(cind);
    cents(rind,:,2) = [cx, cy];
    
    npre = floor(cents(rind,:,2)-cents(rind,:,1));
    ckern = shiftmat(ks(:,:,rind,1),npre(1),npre(2));
    
    ks(:,:,rind,2) = ckern;
    
%     figure(1);clf
%     subplot(1,2,1)
%     showkernel(ks(:,:,rind,1),cents(rind,:,1))
%     subplot(1,2,2)
%     showkernel(ks(:,:,rind,2),cents(rind,:,2))
%     keyboard
end

figure(1);clf
% alsubplot(3*nfigrow,length(fnames)/nfigrow,1,1);
alsubplot(nfigrow,length(fnames)/nfigrow,1,1);
for i = 1:length(fnames)
    cim = im2double(rgb2gray(imread([dname,fnames{i}])));

    [fy,fx] = ind2sub([nfigrow, length(fnames)],i);
    
    alsubplot(fy,fx)
%     imshow(cim);

    rim = circshift(cim,[0 round(size(cim,2)*90/360)]);
%     alsubplot(fy+1,fx)
%     imshow(rim);

    xdiff = (im_size(2)-newsz(2))/2;

    tcim = cim(:,xdiff+1:end-xdiff);
    trim = rim(:,xdiff+1:end-xdiff);
    diffs = NaN(1,size(ks,4));
    for j = 1:size(ks,4)
        diffs(j) = getRMSdiff(getacts(tcim,ks(:,:,:,j)),getacts(trim,ks(:,:,:,j)));
    end

%     alsubplot(fy+2,fx)
    bar(diffs,'w')
    title(fnames{i})
%     set(gca,'XTickLabel',length(kerns))
%     xlabel('number of RFs')
%     ylabel('RMS difference')
    ylim([0 .25])
end