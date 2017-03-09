function drosonav
    wn = 10;
    kwid = 135*2;

    load('../matant/code/idfridf/views/views001.mat');
%     load(sprintf('../../ant/worlds/star%03d.mat',wn),'ims','xg','yg');
%     mx = 2;
%     sel = abs(xg)<=mx/sqrt(2) & abs(yg)<=mx/sqrt(2);
    ims = views(:,:,:);
    xg = mx(:);
    yg = my(:);
%     ims = ims(:,:,sel);
%     xg = xg(sel);
%     yg = yg(sel);

    load('vf_kernels.mat');
    
    main(kwid,ims,xg,yg,vf_avkernels_r2,vf_avkernels_r4,thresh);
end

function main(kwid,ims,xg,yg,vf_avkernels_r2,vf_avkernels_r4,thresh)
    ind = randperm(28);
    acts_r2 = kernelacts(vf_avkernels_r2(ind(1:14)));
    acts_r4 = kernelacts(vf_avkernels_r4);
    
    [~,ref] = min(hypot(xg,yg));
    idf_r2 = myidf(acts_r2);
    idf_r4 = myidf(acts_r4);
    idf_both = myidf([acts_r2,acts_r4]);
    
    showfig(1,idf_r2);
    title('R2 neurons')
    showfig(2,idf_r4);
    title('R4d neurons')
    showfig(3,idf_both);
    title('R2 & R4d neurons')

    function acts=kernelacts(kernels)
        kernels = {kernels.k};
        kernels = kernels(~cellfun(@isempty,kernels));
        acts = NaN(length(xg),length(kernels));
        for i = 1:length(kernels)
            ckern = kernels{i};
            ksz = [size(ims,1) round(size(ims,2)*kwid/360)];
            ckern = resizekernel(ckern,ksz,thresh);
            pad = (size(ims,2)-ksz(2))/2;
            ckern = [zeros(ksz(1),floor(pad)),ckern,zeros(ksz(1),ceil(pad))];
            acts(:,i) = shiftdim(sum(sum(bsxfun(@times,ims,ckern))));
        end
    end

    function idf=myidf(acts)
        idf=sqrt(mean(bsxfun(@minus,acts,acts(ref,:)).^2,2));
    end

    function showfig(fig,idf)
        figure(fig);clf
        [idfim,mx,my] = makeim(xg,yg,idf);
        surf(mx,my,idfim);
        xlabel('x (m)');
        ylabel('y (m)');
        zlabel('Difference in activation');
    end
end