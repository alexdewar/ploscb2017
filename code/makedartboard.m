imsz = [636 2160];
nring = 8.9;
thoff = pi/20;
nseg = 3;
expconst = 2.5;

dartboardim = 255.*ones(imsz(1),imsz(2)/2,'uint8');

[yi,xi] = ndgrid(1:imsz(1),1:(imsz(2)/2));
[th,r] = cart2pol(xi-1,yi-1);

whr = exp(expconst.*linspace(0,1,nring));
whr = [max(r(:)).*(whr-min(whr))./range(whr) inf];
whth = [0 thoff+linspace(0,pi/2,nseg+1)];

for i = 2:length(whr)
    for j = (1+mod(i,2)):2:length(whth)-1
        dartboardim(r > whr(i-1) & r <= whr(i) & th > whth(j) & th <= whth(j+1)) = 0;
    end
end

dartboardim = [dartboardim fliplr(dartboardim)];

figure(1);clf
imshow(dartboardim);