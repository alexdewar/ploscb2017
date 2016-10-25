clear
close all

majoraxis = 1;
minoraxis = 0.5;
im_size = [100 100];
nvar = 10;
freqmax = 5;

thoff = pi/5;

scale = linspace(0,1,nvar);
rat = linspace(0,1,nvar);
freq = round(linspace(1,freqmax,nvar));

% refim = ~ellblob(0,2,0,majoraxis,minoraxis,im_size);
% refrp = regionprops(refim,'MajorAxisLength','MinorAxisLength');
% refmja = refrp.MajorAxisLength;
% refmna = refrp.MinorAxisLength;

[mja,mna] = deal(NaN(nvar));

for i = 1:nvar
    for j = 1:nvar
        im = ~ellblob(scale(i),1,freq(j),0,majoraxis,minoraxis,thoff,im_size);
        rp = regionprops(im,'MajorAxisLength','MinorAxisLength');

        figure(10);clf
        imshow(im);
%         colormap gray
        drawnow
        pause(0.5)
        
        mja(i,j) = rp.MajorAxisLength;
        mna(i,j) = rp.MinorAxisLength;
    end
end

figure(1);clf
subplot(1,3,1)
% [amp,I] = sort(amp);
% mna = mna(I);
% mja = mja(I);
% freq = freq(I);
contourf(scale,freq,mna-mna(1,1))
colorbar

subplot(1,3,2);
contourf(scale,freq,mja-mja(1,1)) %,0,refmja,'r+');
colorbar

subplot(1,3,3);
contourf(scale,freq,mna./mja)
colorbar

% title(num2str(refmja./refmna));