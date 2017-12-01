close all
clear all

%% constants
dosave = true;
fnprefix = '~/nsd/arenavid';
imfn = [mfiledir '/../drosodata/antoinestim/touse/09_3_34_162_1_+232_triangles.png'];
nth = 1000;
degpersec = 100; %deg/s
drumht = 1;
drumr = 0.5;
fps = 30;
vlen = 30;
fov = [120 270];

%% compute all heads
radpersec = degpersec*pi/180;
radperfr = radpersec/fps;

heads = mod(cumsum([0,randn(1,199)*pi/100]),2*pi);
allheads = heads(1);
for i = 1:length(heads)-1
    allheads = [allheads,heads(i):sign(heads(i+1)-heads(i))*radperfr:heads(i+1)];
    if length(allheads)/fps >= vlen
        break;
    end
end

%% image
im = im2double(rgb2gray(imread(imfn)));
imsz = size(im);
% figure(1);clf
% imshow(im)

[X,Y] = bw2polygon(im < 0.5);
% figure(1);clf
% alfill(X,Y,'k')
pz = Y*drumht/imsz(1);

th = linspace(0,2*pi,nth);
[cx,cy] = pol2cart(th,drumr);

xoff = round((imsz(2)-fov(2))/2);
imind = xoff+1:imsz(2)-xoff;
imy = linspace(-fov(1)/2,fov(1)/2,imsz(1));
imx = linspace(-fov(2)/2,fov(2)/2,imsz(2));
% im = im;
imyax = -fov(1)/2:30:fov(1)/2;
imxax = -fov(2)/2:45:fov(2)/2;
% [imx,imy] = ndgrid(imx,imy);

%% kernels
load('vf_kernels.mat','vf_avkernels_r2','vf_avkernels_r4');
kerns = [vf_avkernels_r2,vf_avkernels_r4];
lefts = cell2mat({kerns.isleft});
kerns = [kerns(lefts),kerns(~lefts)];
rkerns = resizekernel(kerns,[imsz(1),fov(2)],.25);
% r2ind = 1:28;
% r4ind = 28+(1:14);

%% video loop
if dosave
    vcnt = 1;
    while true
        vfn = sprintf('%s (%04d).avi',fnprefix,vcnt);
        if ~exist(vfn,'file')
            break;
        end
        vcnt = vcnt+1;
    end
    vwrite = VideoWriter(vfn);
    vwrite.FrameRate = fps;
    vwrite.open;
end
for i = 1:length(allheads)
%     tic
    [px,py] = pol2cart(allheads(i)+pi+2*pi*(X-1)/imsz(2),drumr);

    figure(1)
    subplot(3,1,1)
    cim = circshift(im,[0 round(allheads(i)*imsz(2)/(2*pi))]);
    crim = cim(:,imind);
    image(imx,imy,255*(crim>=0.5))
    colormap gray
    set(gca,'XTick',imxax,'YTick',imyax)
    title('"Fly''s view"');
    axis equal tight on

    subplot(3,1,[2 3])
    hold off
    plot3(cx,cy,zeros(size(cx)),'k',cx,cy,drumht*ones(size(cx)),'k')
    hold on
    alfill(px,py,pz,'k')
    axis equal off
    set(gcf,'Color','w')

    text(-1.5*drumr,0,drumht*1.5,'Laser','Color','r','HorizontalAlignment','Center','VerticalAlignment','Bottom')
    if mod(allheads(i)-pi/8,pi) < pi/4
        line([-1.5*drumr 0],[0 0],drumht*[1.5 0.5],'Color','r','LineWidth',2)
    end

    [lx,ly] = pol2cart(pi/4+(0:pi/2:1.5*pi)+allheads(i),drumr);
    line([lx;lx],[ly;ly],drumht*[0 0 0 0;1 1 1 1],'Color','k')
    
    % tether
    line([0 0],[0 0],drumht*[0.5 1.5],'Color','k','LineWidth',2)
    text(0,0,drumht*1.5,'Tether','HorizontalAlignment','Center')
    
    % fly
    plot3(0,0,drumht/2,'bo','MarkerSize',10,'LineWidth',3)
    
%     subplot(2,2,3)
%     acts_r2 = getacts(crim,rkerns);
% %     acts = ones(1,28);
% %     bar(acts)
% %     ylim([-1 1])
%     imagesc(reshape((acts_r2+1)/2,21,2)')
%     caxis([0 1])
%     title('Ring neuron activity')
%     set(gca,'XTick',[],'YTickLabel',{'','Left','','Right'});
%     colormap hot
%     axis tight
    
%     subplot(2,2,5)
%     acts_r2 = getacts(crim,rkerns(:,:,r4ind));
% %     acts = ones(1,28);
% %     bar(acts)
% %     ylim([-1 1])
%     imagesc(reshape((acts_r2+1)/2,7,2)')
%     caxis([0 1])
%     title('R4 neuron activation')
%     set(gca,'XTick',[],'YTickLabel',{'Left','Right'});
%     colormap hot
%     axis equal tight
    
    if dosave
        vwrite.writeVideo(getframe(gcf));
    else
        drawnow
    end

%     pause(1/fps-toc)
end

if dosave
    vwrite.close;
end