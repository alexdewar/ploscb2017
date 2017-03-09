stim = true(112,252);
wd = 20;
stim(:,round(size(stim,2)/5+(-wd/2:wd/2))) = false;
save('dummystim','stim');

figure(1);clf
imshow(stim)