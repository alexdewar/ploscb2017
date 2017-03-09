function applykernels(imfn,ktypes,gns,ksz,isleft,domirror)
if nargin < 6
    domirror = false;
end
if nargin < 5
    isleft = true;
end
if nargin < 4
    ksz = (1:2:7)/2;
end
if nargin < 3
    gns = [];
end
usegabor=false;
if nargin < 2
    ktypes = {'r2','r4'};
else
    fprintf('file: %s\ncell: %s\nglomerulus: %d\nkernel size: %d\nleft: %d\nmirror: %d\n\n', ...
        imfn,ktypes,gns,ksz,isleft,domirror);
    if strcmp(ktypes,'none')
        usegabor=true;
        its = 1;
    else
        ktypes = {ktypes};
        its = numel(ktypes);
    end
end

if usegabor
    load('gb_kernels.mat');
else
    load('vf_kernels.mat');
end

im = loadim(imfn);
if size(im,3)==3
    im = rgb2gray(im);
end
im = im2double(im);

if domirror
    im = [flipud(im);im];
%     figure(1);clf
%     imshow(im);
end

if nargin==1
    startprogbar(1,2*numel(ksz)*sum(~cellfun(@isempty,{vf_avkernels_r2.k,vf_avkernels_r4.k})));
end
for i = 1:its
    if usegabor
        kernels = gb_kernels;
        whkern = gns;
    else
        ktype = ktypes{i};
        kstruct = eval(['vf_avkernels_' ktype]);
        kstruct = kstruct(cell2mat({kstruct.isleft})==isleft);
        if ~isempty(gns)
            kstruct = kstruct(cell2mat({kstruct.glomnum})==gns);
        end
        kernels = {kstruct.k};
        whkern = 1:numel(kernels);
    end
    
%     im = imresize(im,2);
    for j = whkern
        if ~isempty(kernels{j})
            for k = 1:numel(ksz)
                if usegabor
                    outf = sprintf('./convims/%s_gabor%d_×%.1f.jpg',trimext(imfn),j,ksz(k));
                else
                    outf = sprintf('./convims/%s_%s_g%02d_×%.1f_l%d.jpg',trimext(imfn),ktype,kstruct(j).glomnum,ksz(k),isleft);
                end
                if exist(outf,'file')
                    fprintf('%s already exists!\n',outf);
                else
                    ckern = resizekernel(kernels{j},ksz(k),thresh);
                    cim = conv2(im,ckern,'valid');
                    cim = (1+cim./max(abs(cim(:))))/2;
                    szoff = floor((size(im)-size(cim))/2);
                    outim = NaN(size(im));
                    outim(szoff(1)+(1:size(cim,1)),szoff(2)+(1:size(cim,2))) = cim;
                    disp(outf);
                    imwrite(im2uint8(outim),neuroncolormap,outf);
                    if nargin==1 && progbar
                        return;
                    end
                end
            end
        end
    end
end

% for i = 1:numel(kernels)
%     figure(1);clf
%     imshow(cims(:,:,i));
%     pause(1)
% end