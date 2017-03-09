function patt_bar(dosave)
    if ~nargin
        dosave = false;
    end

    patt1 = rgb2gray(im2double(imread(sprintf('%s/../../data/antoinestim/touse/09_3_34_162_1_+232_triangles.png',mfiledir))));
    patt2 = rgb2gray(im2double(imread(sprintf('%s/../../data/antoinestim/touse/09_4_00_162_1_-107_triangles_com.png',mfiledir))));
    
    load('vf_kernels');
    rkerns = resizekernel(vf_avkernels_r2,[120 270],.25);
    
    figure(1);clf
    xoff = 45;

    acts1 = getacts(patt1(:,xoff+(1:270)),rkerns);
    acts2 = getacts(patt2(:,xoff+(1:270)),rkerns);
    
    h=bar([acts1, acts2]);
    set(h(1),'FaceColor',[0 0 .9]);
    set(h(2),'FaceColor',[.9 0 0]);
    ylim([-1 1]);
    xtl = num2str(([1:14,1:14])','%d');
    set(gca,'FontSize',8);
    set(gca,'YTick',-1:0.5:1,'XTick',1:length(acts1),'XTickLabel',xtl);
    xlim([0 length(acts1)+1])
    
    if dosave
        savefig('patt_bar',[20 4]);
    end
    
%     dump2base(true);
end

% function pbfig(im,rkerns)
%     xoff = 45;
% 
%     im = im2double(im);
%     acts1 = getacts(im(:,xoff+(1:270)),rkerns);
%     imr = circshift(im,[0 90]);
%     acts2 = getacts(imr(:,xoff+(1:270)),rkerns);
%     
%     h=bar([acts1, acts2]);
%     set(h(1),'FaceColor',[0 0 .9]);
%     set(h(2),'FaceColor',[.9 0 0]);
%     ylim([-1 1]);
%     xtl = num2str(([1:14,1:14])','%d');
%     set(gca,'FontSize',8);
%     set(gca,'YTick',-1:0.5:1,'XTick',1:length(acts1),'XTickLabel',xtl);
% end