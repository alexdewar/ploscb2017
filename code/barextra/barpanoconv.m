function barpanoconv(dosave)
if ~nargin
    dosave=false;
end

thresh = 0.25;
imfn = 'bar_vfat.png';
im = im2double(loadim([mfiledir '/' imfn]));
if ndims(im)==3
    im = rgb2gray(im);
end

krng = [-60 60; -135 135]';
irng = [-60 60; -180 180]';

% im = im2double(rgb2gray(imread(imfn)));
load('vf_kernels.mat');
kerns = vf_avkernels_r4;

datafn = sprintf('%s/bar_%s.mat',mfiledir,imfn(1:end-4));
if exist(datafn,'file')
    load(datafn);
else
    if any(size(im) < size(kerns(1).k))
        im = imresize(im,size(kerns(1).k));
        rkerns = cell2mat(shiftdim({kerns.k},-1));
    else
        rkerns = resizekernel(kerns,[size(im,1), size(im,2)*270/360],.25);
    end
    % figure(1);clf
    % imagesc(sign(kernel));
    % colormap(neuroncolormap);
    midi = floor(size(rkerns,2)/2);
    oldimwd = size(im,2);
    vals = NaN(length(kerns),oldimwd);
    % im2 = [im(:,end-midi+1:end) im im(:,1:size(kernel,2)-midi)];
    % figure(1);clf
    % imshow(im)
    % return
    xoff = (size(im,2)-size(rkerns,2))/2;
    startprogbar(10,length(vals))
    for i = 1:length(vals)
        cim = circshift(im,[0 1-i-midi]);
        vals(:,i) = getacts(cim(:,xoff+1:end-xoff),rkerns);
        if progbar
            return;
        end
    %     figure(2);clf
    %     imshow(im2(:,i-1+(1:size(kernel,2))))
    end

    save(datafn,'vals');
end

%%
ths = linspace(irng(1,2),irng(2,2),length(vals));
figure(1);clf
% subplot(1,2,1)
lefts = cell2mat({kerns.isleft});
plot(ths,mean(vals(lefts,:)),'b',ths,mean(vals(~lefts,:)),'b--')
% plot(ths,mean(vals))
axis tight
ylim(.5*[-1 1])
set(gca,'XTick',-180:180:180);

% subplot(1,2,2)
% plot(ths,mean(vals(~lefts,:)),'k--');
% format_ticks;
% ylabel('Activation')
% axis tight

if dosave
    savefig(['bar_' imfn(1:end-4)],[20 10])
end