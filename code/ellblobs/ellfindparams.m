function ellfindparams(starti,nvar,nwave,freqmax,job)
majoraxis = 1;
minoraxis = 0.5;
im_size = [100 100];
% nvar = 10;
% freqmax = 5;

thoff = rand(nvar,1)*2*pi;

scale = rand(nvar,1);
amp = rand(nvar,nwave);
phi = 2*pi*rand(nvar,nwave);
freq = randi(freqmax,nvar,nwave);
% refim = ~ellblob(0,2,0,majoraxis,minoraxis,im_size);
% refrp = regionprops(refim,'MajorAxisLength','MinorAxisLength');
% refmja = refrp.MajorAxisLength;
% refmna = refrp.MinorAxisLength;

[mja,mna] = deal(NaN(nvar,1));

startprogbar(10,nvar);
for i = 1:nvar
    im = ~ellblob(scale(i),amp(i,:),freq(i,:),phi(i,:),majoraxis,minoraxis,thoff(i),im_size);
    rp = regionprops(im,'MajorAxisLength','MinorAxisLength');

    mja(i) = rp.MajorAxisLength;
    mna(i) = rp.MinorAxisLength;
    
    if progbar
        return;
    end
end

fname = sprintf('%s/ellfp%05d_%05d.mat',mfiledir,job,starti);
fprintf('Saving to %s...\n',fname);
save(fname,'majoraxis','minoraxis','im_size','thoff','scale','amp','phi','freq','mna','mja','freqmax','nvar','nwave','job','starti');