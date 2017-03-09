load('vf_kernels.mat');

r2r4 = {'r2';'r4'};

cols = [neuroncolormap(1,:); 1 1 1; neuroncolormap(end,:)];
RL = 'RL';

for i = 1:numel(r2r4)
    kernels = eval(['vf_avkernels_' r2r4{i}]);
    
    for j = 1:numel(kernels)
        ckern = 2+sign(kernels(j).k);
        im = NaN([size(ckern),3]);
        for k = 1:size(im,1)
            im(k,:,:) = cols(ckern(k,:),:);
        end
%         figure(1);clf
%         imshow(ckern);
%         colormap(neuroncolormap);
%         keyboard
        imwrite(im,sprintf('./av_RF_pics/%s_g%02d_%s.png',r2r4{i},kernels(j).glomnum,RL(kernels(j).isleft+1)));
    end
end