function barpanoconv_varywid(dosave)
if ~nargin
    dosave=false;
end

thresh = 0.25;

krng = [-60 60; -135 135]';
irng = [-60 60; -180 180]';

imsz = [range(irng(:,1)), range(irng(:,2))];
ksz = [range(krng(:,1)), range(krng(:,2))];

wids = 30:100:330;

% im = im2double(rgb2gray(imread(imfn)));
load('vf_kernels.mat','vf_avkernels_*');
kerns = [vf_avkernels_r2,vf_avkernels_r4];
rkerns_r2 = resizekernel(vf_avkernels_r2,ksz,thresh);
r2_l = cell2mat({vf_avkernels_r2.isleft});
rkerns_r4 = resizekernel(vf_avkernels_r4,ksz,thresh);
r4_l = cell2mat({vf_avkernels_r4.isleft});

datafn = sprintf('%s/vw_bar.mat',mfiledir);
if exist(datafn,'file')
    load(datafn);
else
    imx = linspace(irng(1,2),irng(2,2),range(irng(:,2)));
    
    [vals_r2_l,vals_r2_r,vals_r4_l,vals_r4_r] = deal(NaN(imsz(2),length(wids)));
    for i = 1:length(wids)
        cim = repmat(abs(imx) > wids(i)/2,imsz(1),1);
        
        imwrite(cim,sprintf('%s/vw_bar%d.png',mfiledir,wids(i)))
        
        fprintf('width: %d deg\n',wids(i))
        disp('R2 L')
        vals_r2_l(:,i) = panoconv_above0(cim,rkerns_r2(:,:,r2_l),ksz(2));
        disp('R2 R')
        vals_r2_r(:,i) = panoconv_above0(cim,rkerns_r2(:,:,~r2_l),ksz(2));
        disp('R4 L')
        vals_r4_l(:,i) = panoconv_above0(cim,rkerns_r4(:,:,r4_l),ksz(2));
        disp('R4 R')
        [vals_r4_r(:,i),ths] = panoconv_above0(cim,rkerns_r4(:,:,~r4_l),ksz(2));
        disp('.')
    end

    save(datafn,'vals*','ths');
end

% endhash

figure(1);clf
for i = 1:size(vals_r2_l,2)
    subplot(1,size(vals_r2_l,2),i)
    plot(ths,vals_r2_l(:,i),ths,vals_r2_r(:,i),'b--', ...
         ths,vals_r4_l(:,i),'g',ths,vals_r4_r(:,i),'g--')
    
%     plot(ths,vals_r2_l(:,i))
%     subplot(4,size(vals_r2_l,2),2)
%     plot(ths,vals_r2_r(:,i),'b--')
%     subplot(4,size(vals_r2_l,2),3)
%     plot(ths,vals_r4_l(:,i),'y')
%     subplot(4,size(vals_r2_l,2),4)
%     plot(ths,vals_r4_r(:,i),'y--')
    
    axis tight
%    ylim(0.5*[-1 1])
    
    if i==1
        set(gca,'FontSize',8,'XTick',-180:180:180); %,'YTick',-0.5:0.5:0.5);
    else
        set(gca,'XTick',[],'YTick',[]);
    end
%     format_ticks(gca);
end
% plot(ths,mean(vals(lefts,:)),'k',ths,mean(vals(~lefts,:)),'k--')


% subplot(1,2,2)
% plot(ths,mean(vals(~lefts,:)),'k--');

% ylabel('Activation')
% axis tight

if dosave
    savefig('vw_bar',[20 10])
end
