clear

load('/home/a/ad/ad374/git/drosopaper/code/drosodata/ANNs/new_ellblob_elaz_job0000000_ind00001.mat')

rim_px_n = prod(pgeb.rim_size,2);
rim_whpx = cumsum([0; rim_px_n]);

allones = NaN(size(x_im,1),3);
for i = 1:3
    allones(:,i) = all(x_im(:,rim_whpx(i)+1:rim_whpx(i+1))==1,2);
end
