clear

patt = rgb2gray(im2double(imread(sprintf('%s/../drosodata/antoinestim/touse/09_4_00_162_1_-107_triangles_com.png',mfiledir))));
load('rx_idf_kerns2.mat','rkernnum','rkerns');
kr1 = rkerns(:,:,rkernnum==1);
kr2 = rkerns(:,:,rkernnum==2);

xoff = (size(patt,2)-size(rkerns,2))/2;

patt1 = patt(:,xoff+(1:size(rkerns,2)));
acts1_r1 = getacts(patt1,kr1);
acts1_r2 = getacts(patt1,kr2);

patt2 = circshift(patt,[0 90]);
patt2 = patt2(:,xoff+(1:size(rkerns,2)));
acts2_r1 = getacts(patt2,kr1);
acts2_r2 = getacts(patt2,kr2);

dr1 = getRMSdiff(acts1_r1,acts2_r1);
dr2 = getRMSdiff(acts1_r2,acts2_r2);

figure(1);clf
bar([dr2 dr1])