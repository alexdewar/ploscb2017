function vf_rotavav

progbaron = false;
ref_kern_i = 1;
nrots = 359;
thresh = 0.25;
minarea = 20;

load('vf_tkernels.mat');
tkerns = tkerns_r2;
tkerns = tkerns(cell2mat({tkerns.isleft}));

szs = cell2mat(cellfun(@size,{tkerns.k}','UniformOutput',false));
im_size = 2*max(max(szs))*[1 1];
tref_kern = sign(tkerns(ref_kern_i).k);
avkernk = zeros(im_size);
avkernk(1:size(tref_kern,1),1:size(tref_kern,2)) = tref_kern;
avkernk = shiftmat(avkernk, im_size(2)/2 - tkerns(ref_kern_i).cent(1), ...
                            im_size(1)/2 - tkerns(ref_kern_i).cent(2));
inh = avkernk==-1;

% figure(1);clf
% showkernel(ref_kern);

rots = linspace(0,360,nrots+1);
rots = rots(1:end-1);
if progbaron
    startprogbar(10,(length(tkerns)-1)*nrots);
end
for i = 3 %:length(tkerns)
    ctkern = sign(tkerns(i).k);
    ckern = zeros(im_size);
    ckern(1:size(ctkern,1),1:size(ctkern,2)) = ctkern;
    ckern = shiftmat(ckern, im_size(2)/2 - tkerns(i).cent(1), ...
                            im_size(1)/2 - tkerns(i).cent(2));
	cinh = ckern==-1;
                        
% 	  figure(2);clf
%     imshow(ckern)
%     keyboard

    bestol = 0;
    bestrot = NaN;
    for j = 1:nrots
        rinh = imrotate(cinh,rots(j),'nearest','crop');
        
%         figure(3);clf
%         imshow(rinh);
%         drawnow;
%         pause(0.5);

        ol = sum(sum(rinh & inh));
        if ol > bestol
            bestrot = j;
            bestol = ol;
        end
        
        if progbaron && progbar
            return;
        end
    end
    
    rcex = imrotate(ckern==1, rots(bestrot),'nearest','crop');
    rcin = imrotate(ckern==-1,rots(bestrot),'nearest','crop');
    
%     rctkern = rcex - rcin;
%     figure(5);clf
%     showkernel(rctkern);
%     keyboard
    
    avkernk = avkernk + rcex - rcin;
end

% avkernk = avkernk/length(tkerns);
%avkernk = (abs(avkernk)>=thresh) .* sign(avkernk);

try
    figure(10);clf
    showkernel(avkernk);
catch ex
    disp(ex)
end

pos = avkernk==1;
neg = avkernk==-1;

bwl = bwlabeln(pos);
bwlm = bwlabeln(neg);
emp = bwl==0;
bwl(emp) = bwlm(emp)+max(bwl(:));

rp = regionprops(bwl,'Area');
areas = cell2mat({rp.Area});

for i = find(areas<minarea)
    avkernk(bwl==i) = 0;
end

avkernk = bwtrim(avkernk);

pos = avkernk==1;
neg = avkernk==-1;
avkernk(pos) = 1/sum(pos(:));
avkernk(neg) = -1/sum(neg(:));

bwlp = bwlabeln(pos);
prp = regionprops(bwlp,'Area','Centroid');
pareas = cell2mat({prp.Area});

[~,whbig] = max(pareas)

avrotkern = struct('k',avkernk,'cent',prp(whbig).Centroid);

vals = unique(bwlp(:));
vals = vals(vals>0);
figure(1);clf
subplot(1,length(vals)+1,1)
showkernel(avrotkern)
for i = 1:length(vals)
    subplot(1,length(vals)+1,i+1)
    imshow(bwlp==vals(i))
end
keyboard

save(fullfile(mfiledir,'vf_rotavav.mat'),'avrotkern');