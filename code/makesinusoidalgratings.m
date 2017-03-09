freqs = [1 2 4 8 16 32];
sz = 500;

for i = 1:numel(freqs)
    fmax = freqs(i)*2*pi;
    
    vals = (sin(fmax*(0:sz-1)/(sz-1))+1)/2;
    imvert = repmat(vals,sz,1);
    save(sprintf('singrid_v_f%02d_sz%d',freqs(i),sz),'imvert');
    
    imhorz = repmat(vals',1,sz);
    save(sprintf('singrid_h_f%02d_sz%d',freqs(i),sz),'imhorz');
    
    imminus45 = NaN(sz);
    for j = 1:sz
        imminus45(j,:) = circshift(vals,[0 j]);
    end
%     figure(1);clf
%     imshow(imminus45);
%     keyboard
    save(sprintf('singrid_-45_f%02d_sz%d',freqs(i),sz),'imminus45');
    
    implus45 = fliplr(imminus45);
    save(sprintf('singrid_+45_f%02d_sz%d',freqs(i),sz),'implus45');
    
%     figure(1);clf
%     imshow(imminus45);
%     keyboard
end
