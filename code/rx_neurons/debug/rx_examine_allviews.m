clear
close all

warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');

fname = 'rx_views_nest1_drum1.000000e-03_z.mat';

nrow = 20;
ncol = 5;

perfig = nrow*ncol;

load([mfiledir '/../../../data/rx_neurons/views/' fname]);

for i = 1:size(views,3)
    crow = 1+mod(i-1,nrow);
    ccol = 1+mod(floor((i-1)./nrow),ncol);
    if crow==1 && ccol==1
        h=figure(1+floor((i-1)./perfig));
        clf
        
        pause(0.01);
        frame_h = get(h,'JavaFrame');
        set(frame_h,'Maximized',1);
    end
    alsubplot(nrow,ncol,crow,ccol)
    imshow(lr_views(:,:,i))
    if ccol==1
        ylabel(crow)
    end
    if crow==nrow
        xlabel(i-nrow)
    end
end