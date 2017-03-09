% function changeweights
clear

doprogbar = true;
dosave = false;

pattthrng = 30;

frows = 2;
fcols = 4;

imsz = [120 360];
fov = [120 270];
ksz = fov;
figsz = [20 20];
dname = [mfiledir '/../dispatchrecap/antoinestim/touse'];

load('vf_kernels.mat','vf_avkernels_r2');
kerns = vf_avkernels_r2;
rkerns = resizekernel(kerns,ksz,0.25);

d = dir(fullfile(dname,'*.png'));
if doprogbar
    startprogbar(1,length(d))
end
fi = NaN(size(d));
for i = 1:length(fi)
    fi(i) = (d(i).name(6)~='0') + 2*(d(i).name(7)~='0');
end
[~,I] = sort(fi);
d = d(I);
pcnt = 1;
for i = 1:length(d)
    whrow = ceil(i/fcols);
    whcol = 1+mod(i-1,fcols);
    whfig = ceil(whrow/frows);
    whrow = 1+mod(whrow-1,frows);
    if whrow==1 && whcol==1
        figure(whfig);clf
        alsubplot(frows*2,fcols,1,1);
        title(num2str(whfig))
        hold on
    end
    
    im = rgb2gray(im2double(imread(fullfile(dname,d(i).name))));
    [acts,angs] = panoconv_all(im,rkerns,fov(2));
    Xtrain = [ones(size(acts,2),1), acts'];
%     patt1 = find(abs(angs)<=pattthrng/2);
%     patt2 = find(abs(angs-90)<=pattthrng/2);
%     Xtrain = [acts(:,patt1), acts(:,patt2)]';
%     Xtrain = [ones(size(Xtrain,1),1), Xtrain]; % add bias
    
%     T = [ones(length(patt1),1);-ones(length(patt2),1)];
    T = cosd(2*angs)';
    
    W = pinv(Xtrain)*T;
    
%     mid1 = ceil(length(patt1)/2);
%     mid2 = ceil(length(patt2)/2);
%     Xtest = Xtrain(mid1+[0 mid2-1],:);
    Y = W'*Xtrain';
    
%     [wacts,angs] = wpanoconv(im,rkerns,W,fov(2));

%     keyboard
%     weights = 1./act0;
%     infs = isinf(weights);
%     weights(infs) = max(weights(~infs));
%     wacts = wpanoconv(im,rkerns,weights,fov(2));
    
    alsubplot(1+2*(whrow-1),whcol);
    imshow(im);
    title(d(i).name([1 2]))
%     
    alsubplot(2+2*(whrow-1),whcol);
%     [~,i90] = min(abs(angs+90));
%     [~,i0] = min(abs(angs));
%     plot(angs,acts,angs,wacts) %, ...
% %         angs([i0,i90]),wacts([i0,i90]),'ro');
    plot(angs,Y,angs,T',angs,sign(Y),'r--')
%     plot(1:length(Y),Y,1:length(Y),T',1:length(Y),sign(Y))
    ylim(1.5*[-1 1])
    xlim([angs(1) angs(end)]);
    if d(i).name(6)=='0'
        spp = '-';
    else
        spp = sprintf('+ (%s)',d(i).name(6));
    end
    if d(i).name(7)=='0'
        dcp = '-';
    else
        dcp = sprintf('+ (%s)',d(i).name(7));
    end
    title(['SPP' spp ' / DCP' dcp])
    set(gca,'XTick',[])
    xlabel(sprintf('e=%f',sum(abs(W))))
%     xlabel(sprintf('diff: %.3f',wacts(i90)-wacts(i0)));
    if doprogbar && progbar
        return
    end
    
    if dosave && (i==length(d) || (whrow==frows && whcol==fcols))
        cnt = savefig(sprintf('patterns_werr_%02d',whfig),figsz);
        pcnt = max(pcnt,cnt);
    end
end

if dosave
    close all
end

% if dosave
%     close all
%     cmd=sprintf('pdftk "%s/figures/%04d_patterns_*" cat output "%s/figures/patterns_%04d.pdf"', ...
%                    mfiledir,pcnt,mfiledir,pcnt);
% 	disp(cmd);
%     system(cmd);
% end