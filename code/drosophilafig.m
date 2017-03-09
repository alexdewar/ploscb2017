function drosophilafig(dosave)
kernelscale = 1/5;
imfiles = { 'cameraman.tif' }; %; 'pears.png'; 'tire.tif'; 'circlesBrightDark.png' };

close all
warning('off','images:imshow:magnificationMustBeFitForDockedFigure');
load('kernels.mat')

% figure(4);
% for i = 1:numel(tkernels)
%     alsubplot(numel(kernels)+1,numel(imfiles)+1,i+1,1);
%     imagesc(tkernels{i});
% %     title(labels{i});
%     colormap(drosophilacolormap);
%     freezeColors;
%     axis equal off
% end

startprogbar(1,numel(imfiles)*numel(tkernels));
% for i = 1:numel(imfiles)
% i = 1;
    im = imread(imfiles{1});
    if ndims(im)==3
        im = rgb2gray(im);
    end
    
%     alsubplot(numel(tkernels)+1,2,1,2);
%     imshow(im);
    bigim = imresize(im2double(im),1/kernelscale);
    bigsz = size(bigim)/4;
    for j = 1:numel(tkernels)
        figure(j);clf
        subplot(1,2,1);
        imagesc(tkernels{j});
        title(labels{j});
        colormap(drosophilacolormap);
        freezeColors;
        axis equal tight

        subplot(1,2,2);
        imshow(im);
        title(labels{j});
        convim = conv2(bigim,tkernels{j});
        imshow(convim);
        colormap gray
        freezeColors;
        hold on
        sz = size(tkernels{j});
        fill(bigsz(2)+sz(2)*[-1 -1 1 1 -1],bigsz(1)+sz(1)*[-1 1 1 -1 -1],'g','LineStyle','none'); %,'FaceAlpha',0.5);
        axis square tight off
        
        progbar;
        
        if dosave
            savefig('drosophila',[4,2]);
        end
    end
% end

% if dosave
%     figure(4);savefig('drosophilaimages',[25 120]);
% end