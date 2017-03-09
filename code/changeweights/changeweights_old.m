% function changeweights
clear

doprogbar = true;
dosave = true;

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
%         title(num2str(whfig))
%         hold on
    end
    
    im = rgb2gray(im2double(imread(fullfile(dname,d(i).name))));
    acts = panoconv_all(im,rkerns,ones(length(kerns),1),fov(2));
%     X = [ones(1,2); act0, act90]';
    actblank = (getacts(rand(ksz),rkerns)+1)/2;
%     actblank = 1-act0;
    X = [act0, actblank]';
    XT = pinv(X);
    W = XT*[0;1];
%     Y = W'*act90;
    [wacts,angs] = wpanoconv(im,rkerns,W,fov(2));

%     keyboard
%     weights = 1./act0;
%     infs = isinf(weights);
%     weights(infs) = max(weights(~infs));
%     wacts = wpanoconv(im,rkerns,weights,fov(2));
    
    alsubplot(1+2*(whrow-1),whcol);
    imshow(im);
    title(d(i).name([1 2]))
    
    alsubplot(2+2*(whrow-1),whcol);
    [~,i90] = min(abs(angs+90));
    [~,i0] = min(abs(angs));
    plot(angs,acts,angs,wacts) %, ...
%         angs([i0,i90]),wacts([i0,i90]),'ro');
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
    xlabel(sprintf('diff: %.3f',wacts(i90)-wacts(i0)));
    if doprogbar && progbar
        return
    end
    
    if dosave && (i==length(d) || (whrow==frows && whcol==fcols))
        cnt = savefig(sprintf('patterns_%02d',whfig),figsz);
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