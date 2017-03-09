
green = zeros(100,100,4);
green(:,:,2) = 1;
green(20:40,20:40,4) = .75;
green(40:80,40:80,4) = .5;

im = imageadd(zeros(100,100),'alpha',0.25,green);

figure(1);clf
imshow(im)