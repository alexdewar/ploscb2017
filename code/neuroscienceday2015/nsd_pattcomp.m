function nsd_pattcomp

imfn1 = '09_3_34_162_1_+232_triangles.png';
imfn2 = '09_4_00_162_1_-107_triangles_com.png';
figsz = [9 5];

load('vf_kernels','vf_avkernels_r2');
rkerns = resizekernel(vf_avkernels_r2,[120 270],0.25);

close all

figure(1);clf
pcacts(imfn1,rkerns);
savefig('nsd_pattcomp1',figsz,'eps')

figure(2);clf
pcacts(imfn2,rkerns);
savefig('nsd_pattcomp2',figsz,'eps')

dump2base(true)

end

function acts=pcacts(imfn,rkerns)
    im = im2double(rgb2gray(imread([mfiledir '/../../data/antoinestim/touse/' imfn])));
    
    acts1=getacts(im(:,45+(1:270)),rkerns);
    rim = circshift(im,[0 90]);
    acts2=getacts(rim(:,45+(1:270)),rkerns);
    
    h=bar(0.5*(1+[acts1, acts2]));
    set(h(1),'FaceColor','w');
    set(h(2),'FaceColor','k');
    ylim([0 .7])
    xlim([0 29])
    
    ylabel('Activation')
    set(gca,'XTick',4:4:28)
    xlabel('R2 filter')
    
    fprintf('===== DIFF: %f =====\n',getRMSdiff(acts1,acts2))
end