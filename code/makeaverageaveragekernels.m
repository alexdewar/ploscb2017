function makeaverageaveragekernels
    load vf_kernels.mat
    avkr2 = doaak(vf_avkernels_r2);
    avkr4 = doaak(vf_avkernels_r4);
    avkr2_left = doaak(vf_avkernels_r2(cell2mat({vf_avkernels_r2.isleft})));
    avkr4_left = doaak(vf_avkernels_r4(cell2mat({vf_avkernels_r4.isleft})));
    
    clear vf*
    save('vf_averageaveragekernels');
end

function avk=doaak(kernels)
    cents = cell2mat({kernels.cent}');
    avcent = mean(cents);
    scents = bsxfun(@minus,avcent,cents);
    avkk = zeros(size(kernels(1).k));
    for i = 1:numel(kernels)
        sk = shiftmat(sign(kernels(i).k),scents(i,1),scents(i,2));
        avkk = avkk + sk(1:size(avkk,1),1:size(avkk,2));
    end
    
    avkk = bwtrim(threshkern(avkk/numel(kernels),0.25));
    rp = regionprops(avkk>0,'Centroid');
%     figure(1);clf
%     showkernel(avkk)
%     keyboard
    
    avk = struct('k',avkk,'cent',rp.Centroid);
end